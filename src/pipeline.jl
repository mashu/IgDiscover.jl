# Pipeline orchestration — runs the full IgDiscover workflow

"""
    init_analysis(dir, db_dir, reads)

Initialize an analysis directory with database, reads, and default config.
"""
function init_analysis(dir::AbstractString, db_dir::AbstractString,
                      reads::AbstractString)
    for subdir in ("database", "reads", "stats")
        mkpath(joinpath(dir, subdir))
    end

    for gene in ("V", "D", "J")
        src = joinpath(db_dir, "$gene.fasta")
        isfile(src) || error("Missing $src")
        cp(src, joinpath(dir, "database", "$gene.fasta"); force = true)
    end

    dst_reads = joinpath(dir, endswith(reads, ".gz") ? "reads.fasta.gz" : "reads.fasta")
    cp(reads, dst_reads; force = true)

    write_default_config(joinpath(dir, "igdiscover.toml"))
    @info "Initialized analysis at $dir"
end

"""
    run_pipeline(dir)

Run the complete IgDiscover pipeline in the given analysis directory.
"""
function run_pipeline(dir::AbstractString)
    cd(dir) do
        config = load_config()
        @info "IgDiscover.jl — $(config.iterations) iteration(s)"
        config.limit > 0 && @info "Read limit active: only first $(config.limit) reads will be used (set limit=0 in igdiscover.toml for full run)"

        reads_path = find_reads()
        if config.barcode_length != 0
            @info "PCR bias correction: grouping by barcode (barcode_length=$(config.barcode_length))"
            reads_path = group_reads(reads_path, "grouped.fasta.gz", config; limit = config.limit)
        end

        for iteration in 1:config.iterations
            iter_dir = @sprintf("iteration-%02d", iteration)
            @info "━━━ Iteration $iteration ━━━"
            run_iteration(config, iter_dir, reads_path, iteration)
        end

        run_final(config, reads_path)
        @info "Pipeline complete"
    end
end

function find_reads()
    for candidate in ("reads/sequences.fasta.gz", "reads.fasta.gz", "reads.fastq.gz")
        isfile(candidate) && return candidate
    end
    error("No reads found (expected reads.fasta.gz or reads/sequences.fasta.gz)")
end

function run_iteration(config::Config, iter_dir::String, reads_path::String, iteration::Int)
    mkpath(joinpath(iter_dir, "database"))
    mkpath(joinpath(iter_dir, "stats"))

    # Set up database
    v_src = if iteration == 1
        joinpath("database", "V.fasta")
    else
        prev = @sprintf("iteration-%02d", iteration - 1)
        joinpath(prev, "new_V_pregermline.fasta")
    end
    cp(v_src, joinpath(iter_dir, "database", "V.fasta"); force = true)
    cp(joinpath("database", "D.fasta"), joinpath(iter_dir, "database", "D.fasta"); force = true)

    j_src = if iteration > 1 && config.j_discovery.propagate && isfile("iteration-01/new_J.fasta")
        "iteration-01/new_J.fasta"
    else
        joinpath("database", "J.fasta")
    end
    cp(j_src, joinpath(iter_dir, "database", "J.fasta"); force = true)

    db_dir = joinpath(iter_dir, "database")

    # IgBLAST
    airr_path = joinpath(iter_dir, "airr.tsv.gz")
    if !isfile(airr_path)
        run_igblast_on_fasta(db_dir, reads_path, airr_path;
            sequence_type = config.sequence_type,
            species = config.species,
            penalty = config.mismatch_penalty,
            threads = Sys.CPU_THREADS,
            limit = config.limit)
    end

    # Augment
    assigned_path = joinpath(iter_dir, "assigned.tsv.gz")
    if !isfile(assigned_path)
        @info "Augmenting..."
        airr_df = read_assignments(airr_path)
        augmented = augment_table(airr_df, db_dir; sequence_type = config.sequence_type)
        write_table_gzipped(assigned_path, augmented)
    end

    # Filter
    filtered_path = joinpath(iter_dir, "filtered.tsv.gz")
    if !isfile(filtered_path)
        filter_table(assigned_path, filtered_path, config.preprocessing_filter;
                    stats_path = joinpath(iter_dir, "stats", "filtered.json"))
    end

    # Discover
    candidates_path = joinpath(iter_dir, "candidates.tab")
    if !isfile(candidates_path)
        table = read_assignments(filtered_path)
        database = read_fasta_dict(joinpath(db_dir, "V.fasta"))
        candidates_df = discover_germline(table, database, config)
        write_table(candidates_path, candidates_df)
    end

    # Germline filters
    is_last = iteration == config.iterations
    run_germline_filters(config, iter_dir, is_last)

    # J discovery (iteration 1 only, simplified — keep existing J)
    if iteration == 1
        j_out = joinpath(iter_dir, "new_J.fasta")
        isfile(j_out) || cp(joinpath(db_dir, "J.fasta"), j_out; force = true)
    end
end

function run_germline_filters(config::Config, iter_dir::String, is_last::Bool)
    candidates_path = joinpath(iter_dir, "candidates.tab")
    wl_paths = String[]
    config.germline_filter.whitelist && isfile("database/V.fasta") && push!(wl_paths, "database/V.fasta")
    isfile("whitelist.fasta") && push!(wl_paths, "whitelist.fasta")
    wl = Whitelist(wl_paths)

    # Always produce both pregermline and germline outputs
    for (prefix, criteria) in [("pre", config.pre_germline_filter), ("", config.germline_filter)]
        fasta_path = joinpath(iter_dir, "new_V_$(prefix)germline.fasta")
        isfile(fasta_path) && continue

        candidates = CSV.read(candidates_path, DataFrame; delim = '\t')
        if nrow(candidates) == 0 || !hasproperty(candidates, :consensus)
            @info "No candidates to filter, writing empty germline outputs"
            write_table(joinpath(iter_dir, "new_V_$(prefix)germline.tab"), candidates)
            write_table(joinpath(iter_dir, "annotated_V_$(prefix)germline.tab"), candidates)
            write_fasta(fasta_path, FastaRecord[])
        else
            filtered, annotated = germline_filter!(candidates, criteria; whitelist = wl)
            write_table(joinpath(iter_dir, "new_V_$(prefix)germline.tab"), filtered)
            write_table(joinpath(iter_dir, "annotated_V_$(prefix)germline.tab"), annotated)
            germline_filter_to_fasta(filtered, fasta_path)
        end
    end
end

function run_final(config::Config, reads_path::String)
    final_dir = "final"
    mkpath(joinpath(final_dir, "database"))
    mkpath(joinpath(final_dir, "stats"))

    last_iter = @sprintf("iteration-%02d", max(config.iterations, 1))

    # Copy final databases (fall back to initial V if discovered set is empty)
    v_src = if config.iterations > 0 && isfile(joinpath(last_iter, "new_V_germline.fasta"))
        cand = joinpath(last_iter, "new_V_germline.fasta")
        isempty(read_fasta(cand)) ? joinpath("database", "V.fasta") : cand
    else
        joinpath("database", "V.fasta")
    end
    isfile(v_src) || error("V database not found at $v_src")
    cp(v_src, joinpath(final_dir, "database", "V.fasta"); force = true)
    cp(joinpath("database", "D.fasta"), joinpath(final_dir, "database", "D.fasta"); force = true)

    j_src = if config.j_discovery.propagate && isfile("iteration-01/new_J.fasta")
        "iteration-01/new_J.fasta"
    else
        joinpath("database", "J.fasta")
    end
    cp(j_src, joinpath(final_dir, "database", "J.fasta"); force = true)

    db_dir = joinpath(final_dir, "database")
    airr_path = joinpath(final_dir, "airr.tsv.gz")
    assigned_path = joinpath(final_dir, "assigned.tsv.gz")
    filtered_path = joinpath(final_dir, "filtered.tsv.gz")

    # If we fell back to original V, reuse iteration results to avoid re-running IgBLAST
    if v_src == joinpath("database", "V.fasta") && isfile(joinpath(last_iter, "airr.tsv.gz"))
        @info "Using iteration results for final (no new V genes discovered)"
        cp(joinpath(last_iter, "airr.tsv.gz"), airr_path; force = true)
        cp(joinpath(last_iter, "assigned.tsv.gz"), assigned_path; force = true)
        cp(joinpath(last_iter, "filtered.tsv.gz"), filtered_path; force = true)
        if isfile(joinpath(last_iter, "stats", "filtered.json"))
            cp(joinpath(last_iter, "stats", "filtered.json"), joinpath(final_dir, "stats", "filtered.json"); force = true)
        end
    else
        # Final IgBLAST + augment + filter
        if !isfile(airr_path)
            run_igblast_on_fasta(db_dir, reads_path, airr_path;
                sequence_type = config.sequence_type, species = config.species,
                penalty = config.mismatch_penalty, threads = Sys.CPU_THREADS,
                limit = config.limit)
        end
        if !isfile(assigned_path)
            airr_df = read_assignments(airr_path)
            augmented = augment_table(airr_df, db_dir; sequence_type = config.sequence_type)
            write_table_gzipped(assigned_path, augmented)
        end
        if !isfile(filtered_path)
            filter_table(assigned_path, filtered_path, config.preprocessing_filter;
                        stats_path = joinpath(final_dir, "stats", "filtered.json"))
        end
    end

    @info "Final output in $final_dir/"
end
