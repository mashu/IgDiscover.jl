# Pipeline orchestration — runs the full IgDiscover workflow
#
# Each processing step is a standalone function. File-level caching is
# handled by `cached(path) do ... end` so the pipeline can resume after
# interruption without re-running completed steps.

# ─── File-level caching ───

"""
    cached(f, path) -> path

Execute `f()` only if `path` does not exist. Returns `path`.
Provides a uniform, declarative caching pattern for pipeline steps.
"""
function cached(f::Function, path::AbstractString)
    isfile(path) || f()
    path
end

# ─── Initialization ───

"""
    init_analysis(dir, db_dir, reads)

Initialize an analysis directory with database, reads symlink, and default config.
"""
function init_analysis(dir::AbstractString, db_dir::AbstractString, reads::AbstractString)
    for subdir in ("database", "stats")
        mkpath(joinpath(dir, subdir))
    end
    for gene in ("V", "D", "J")
        src = joinpath(db_dir, "$gene.fasta")
        isfile(src) || error("Missing $src")
        cp(src, joinpath(dir, "database", "$gene.fasta"); force=true)
    end
    dst_reads = joinpath(dir, endswith(reads, ".gz") ? "reads.fasta.gz" : "reads.fasta")
    cp(reads, dst_reads; force=true)
    write_default_config(joinpath(dir, "igdiscover.toml"))
    @info "Initialized analysis at $dir"
end

# ─── Pipeline entry point ───

"""
    run_pipeline(dir)

Run the complete IgDiscover pipeline in the given analysis directory.
"""
function run_pipeline(dir::AbstractString)
    cd(dir) do
        config = load_config()
        @info "IgDiscover.jl — $(config.iterations) iteration(s)"
        config.limit > 0 && @info "Read limit: $(config.limit) (set limit=0 for full run)"

        t_start = time()
        reads_path = find_reads()

        reads_path = preprocess_reads(config, reads_path)

        for iteration in 1:config.iterations
            t_iter = time()
            @info "━━━ Iteration $iteration ━━━"
            run_iteration(config, @sprintf("iteration-%02d", iteration), reads_path, iteration)
            @info @sprintf("Iteration %d completed in %.1f seconds", iteration, time() - t_iter)
        end

        run_final(config, reads_path)
        @info @sprintf("Pipeline complete in %.1f seconds", time() - t_start)
    end
end

function find_reads()
    for path in ("reads/sequences.fasta.gz", "reads.fasta.gz", "reads.fastq.gz")
        isfile(path) && return path
    end
    error("No reads found (expected reads.fasta.gz or reads/sequences.fasta.gz)")
end

function preprocess_reads(config::Config, reads_path::String)
    if config.barcode_length != 0
        @info "PCR bias correction active (barcode_length=$(config.barcode_length))"
        group_reads(reads_path, "grouped.fasta.gz", config; limit=config.limit)
    elseif config.race_g
        trim_reads_race_g(reads_path, "trimmed.fasta.gz"; limit=config.limit)
    else
        reads_path
    end
end

# ─── Iteration database setup ───

function setup_iteration_database(iter_dir::String, config::Config, iteration::Int)
    db_dst = joinpath(iter_dir, "database")

    v_src = iteration == 1 ?
        joinpath("database", "V.fasta") :
        joinpath(@sprintf("iteration-%02d", iteration - 1), "new_V_pregermline.fasta")

    j_src = if iteration > 1 && config.j_discovery.propagate &&
               isfile("iteration-01/new_J.fasta")
        "iteration-01/new_J.fasta"
    else
        joinpath("database", "J.fasta")
    end

    if iteration == 1
        write_sanitized_imgt(v_src, joinpath(db_dst, "V.fasta"))
        write_sanitized_imgt(joinpath("database", "D.fasta"), joinpath(db_dst, "D.fasta"))
        write_sanitized_imgt(j_src, joinpath(db_dst, "J.fasta"))
    else
        cp(v_src, joinpath(db_dst, "V.fasta"); force=true)
        cp(joinpath("database", "D.fasta"), joinpath(db_dst, "D.fasta"); force=true)
        cp(j_src, joinpath(db_dst, "J.fasta"); force=true)
    end

    joinpath(iter_dir, "database")
end

# ─── Individual pipeline steps ───

function step_igblast(config::Config, db_dir::String, reads_path::String, output_path::String)
    cached(output_path) do
        t = time()
        run_igblast_on_fasta(db_dir, reads_path, output_path;
            sequence_type=config.sequence_type, species=config.species,
            penalty=config.mismatch_penalty, threads=Sys.CPU_THREADS,
            limit=config.limit)
        @info @sprintf("IgBLAST: %.1fs", time() - t)
    end
end

function step_augment(config::Config, db_dir::String, airr_path::String, output_path::String)
    cached(output_path) do
        t = time()
        @info "Augmenting..."
        augmented = augment_table(read_assignments(airr_path), db_dir;
            sequence_type=config.sequence_type)
        write_table_gz(output_path, augmented)
        @info @sprintf("Augment: %.1fs", time() - t)
    end
end

function step_filter(config::Config, assigned_path::String, output_path::String, stats_path::String)
    cached(output_path) do
        t = time()
        filter_table(assigned_path, output_path, config.preprocessing_filter;
            stats_path=stats_path)
        @info @sprintf("Filter: %.1fs", time() - t)
    end
end

function step_discover(config::Config, db_dir::String, filtered_path::String, output_path::String)
    cached(output_path) do
        t = time()
        table    = read_assignments(filtered_path)
        database = read_fasta_dict(joinpath(db_dir, "V.fasta"))
        write_table(output_path, discover_germline(table, database, config))
        @info @sprintf("Discovery: %.1fs", time() - t)
    end
end

function step_germline_filters(config::Config, iter_dir::String)
    candidates_path = joinpath(iter_dir, "candidates.tsv")
    wl_paths = String[]
    config.germline_filter.whitelist && push!(wl_paths, "database/V.fasta")
    isfile("whitelist.fasta") && push!(wl_paths, "whitelist.fasta")
    wl = Whitelist(wl_paths)

    for (prefix, criteria) in [("pre", config.pre_germline_filter), ("", config.germline_filter)]
        fasta_path = joinpath(iter_dir, "new_V_$(prefix)germline.fasta")
        isfile(fasta_path) && continue

        candidates = CSV.read(candidates_path, DataFrame; delim='\t')
        if nrow(candidates) == 0 || !hasproperty(candidates, :consensus)
            @info "No candidates; writing empty germline outputs"
            write_table(joinpath(iter_dir, "new_V_$(prefix)germline.tsv"), candidates)
            write_table(joinpath(iter_dir, "annotated_V_$(prefix)germline.tsv"), candidates)
            write_fasta(fasta_path, FastaRecord[])
        else
            filtered, annotated = germline_filter!(candidates, criteria; whitelist=wl)
            filtered = deduplicate_by_consensus(filtered)
            write_table(joinpath(iter_dir, "new_V_$(prefix)germline.tsv"), filtered)
            write_table(joinpath(iter_dir, "annotated_V_$(prefix)germline.tsv"), annotated)
            germline_filter_to_fasta(filtered, fasta_path)
        end
    end
end

function step_rename(config::Config, iter_dir::String)
    config.rename || return
    germline_fasta = joinpath(iter_dir, "new_V_germline.fasta")
    isfile(germline_fasta) && !isempty(read_fasta(germline_fasta)) || return
    renamed_path = joinpath(iter_dir, "new_V_germline_renamed.fasta")
    rename_genes(germline_fasta, joinpath("database", "V.fasta"), renamed_path)
    cp(renamed_path, germline_fasta; force=true)
end

function step_j_discovery(config::Config, db_dir::String, filtered_path::String, iter_dir::String)
    j_out = joinpath(iter_dir, "new_J.fasta")
    cached(j_out) do
        table      = read_assignments(filtered_path)
        j_database = read_fasta_dict(joinpath(db_dir, "J.fasta"))
        discover_j_to_fasta(table, j_database, config, j_out)
        if isempty(read_fasta(j_out))
            @info "J discovery produced no candidates; keeping original J database"
            cp(joinpath(db_dir, "J.fasta"), j_out; force=true)
        end
    end
end

function step_clonotypes(filtered_path::String, output_path::String)
    cached(output_path) do
        isfile(filtered_path) || return
        t = time()
        table = read_assignments(filtered_path)
        nrow(table) == 0 && return
        hasproperty(table, :cdr3) && hasproperty(table, :v_call) && hasproperty(table, :j_call) || return
        ct_df, _ = call_clonotypes(table)
        write_table(output_path, ct_df)
        @info @sprintf("Clonotypes: %d from %d rows in %.1fs", nrow(ct_df), nrow(table), time() - t)
    end
end

# ─── Iteration orchestration ───

function run_iteration(config::Config, iter_dir::String, reads_path::String, iteration::Int)
    mkpath(joinpath(iter_dir, "database"))
    mkpath(joinpath(iter_dir, "stats"))

    db_dir = setup_iteration_database(iter_dir, config, iteration)

    airr_path       = joinpath(iter_dir, "airr.tsv.gz")
    assigned_path   = joinpath(iter_dir, "assigned.tsv.gz")
    filtered_path   = joinpath(iter_dir, "filtered.tsv.gz")
    candidates_path = joinpath(iter_dir, "candidates.tsv")
    stats_path      = joinpath(iter_dir, "stats", "filtered.json")

    step_igblast(config, db_dir, reads_path, airr_path)
    step_augment(config, db_dir, airr_path, assigned_path)
    step_filter(config, assigned_path, filtered_path, stats_path)
    step_discover(config, db_dir, filtered_path, candidates_path)
    step_germline_filters(config, iter_dir)
    step_rename(config, iter_dir)

    iteration == 1 && step_j_discovery(config, db_dir, filtered_path, iter_dir)
end

# ─── Final run ───

function run_final(config::Config, reads_path::String)
    final_dir = "final"
    mkpath(joinpath(final_dir, "database"))
    mkpath(joinpath(final_dir, "stats"))

    last_iter = @sprintf("iteration-%02d", max(config.iterations, 1))

    # Resolve final V database
    germline_fasta = joinpath(last_iter, "new_V_germline.fasta")
    v_src = if isfile(germline_fasta) && !isempty(read_fasta(germline_fasta))
        germline_fasta
    else
        @info "No new V germline sequences discovered; using original database"
        joinpath("database", "V.fasta")
    end

    # Resolve final J database
    j_src = config.j_discovery.propagate && isfile("iteration-01/new_J.fasta") ?
        "iteration-01/new_J.fasta" : joinpath("database", "J.fasta")

    # Copy final databases
    cp(v_src, joinpath(final_dir, "database", "V.fasta"); force=true)
    cp(joinpath("database", "D.fasta"), joinpath(final_dir, "database", "D.fasta"); force=true)
    cp(j_src, joinpath(final_dir, "database", "J.fasta"); force=true)

    db_dir        = joinpath(final_dir, "database")
    airr_path     = joinpath(final_dir, "airr.tsv.gz")
    assigned_path = joinpath(final_dir, "assigned.tsv.gz")
    filtered_path = joinpath(final_dir, "filtered.tsv.gz")
    stats_path    = joinpath(final_dir, "stats", "filtered.json")

    # Reuse iteration results if V database unchanged
    if v_src == joinpath("database", "V.fasta") && isfile(joinpath(last_iter, "airr.tsv.gz"))
        @info "Reusing iteration results for final (V unchanged)"
        for fname in ("airr.tsv.gz", "assigned.tsv.gz", "filtered.tsv.gz")
            cp(joinpath(last_iter, fname), joinpath(final_dir, fname); force=true)
        end
        src_stats = joinpath(last_iter, "stats", "filtered.json")
        isfile(src_stats) && cp(src_stats, stats_path; force=true)
    else
        step_igblast(config, db_dir, reads_path, airr_path)
        step_augment(config, db_dir, airr_path, assigned_path)
        step_filter(config, assigned_path, filtered_path, stats_path)
    end

    @info "Final output in $final_dir/"

    step_clonotypes(filtered_path, joinpath(final_dir, "clonotypes.tsv"))
end
