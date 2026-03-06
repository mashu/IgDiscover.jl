# Pipeline orchestration — runs the full IgDiscover workflow

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

"""
    run_pipeline(dir)

Run the complete IgDiscover pipeline in the given analysis directory.
"""
function run_pipeline(dir::AbstractString)
    cd(dir) do
        config = load_config()
        @info "IgDiscover.jl — $(config.iterations) iteration(s)"
        config.limit > 0 && @info "Read limit: $(config.limit) (set limit=0 for full run)"

        reads_path = find_reads()

        if config.barcode_length != 0
            @info "PCR bias correction active (barcode_length=$(config.barcode_length))"
            reads_path = group_reads(reads_path, "grouped.fasta.gz", config; limit=config.limit)
        end

        for iteration in 1:config.iterations
            @info "━━━ Iteration $iteration ━━━"
            run_iteration(config, @sprintf("iteration-%02d", iteration), reads_path, iteration)
        end

        run_final(config, reads_path)
        @info "Pipeline complete"
    end
end

function find_reads()
    for path in ("reads/sequences.fasta.gz", "reads.fasta.gz", "reads.fastq.gz")
        isfile(path) && return path
    end
    error("No reads found (expected reads.fasta.gz or reads/sequences.fasta.gz)")
end

function run_iteration(config::Config, iter_dir::String, reads_path::String, iteration::Int)
    mkpath(joinpath(iter_dir, "database"))
    mkpath(joinpath(iter_dir, "stats"))

    # Set up per-iteration database
    v_src = iteration == 1 ?
        joinpath("database", "V.fasta") :
        joinpath(@sprintf("iteration-%02d", iteration - 1), "new_V_pregermline.fasta")
    cp(v_src, joinpath(iter_dir, "database", "V.fasta"); force=true)
    cp(joinpath("database", "D.fasta"), joinpath(iter_dir, "database", "D.fasta"); force=true)

    j_src = if iteration > 1 && config.j_discovery.propagate &&
               isfile("iteration-01/new_J.fasta")
        "iteration-01/new_J.fasta"
    else
        joinpath("database", "J.fasta")
    end
    cp(j_src, joinpath(iter_dir, "database", "J.fasta"); force=true)

    db_dir       = joinpath(iter_dir, "database")
    airr_path    = joinpath(iter_dir, "airr.tsv.gz")
    assigned_path = joinpath(iter_dir, "assigned.tsv.gz")
    filtered_path = joinpath(iter_dir, "filtered.tsv.gz")
    candidates_path = joinpath(iter_dir, "candidates.tab")

    # IgBLAST
    if !isfile(airr_path)
        run_igblast_on_fasta(db_dir, reads_path, airr_path;
            sequence_type=config.sequence_type, species=config.species,
            penalty=config.mismatch_penalty, threads=Sys.CPU_THREADS,
            limit=config.limit)
    end

    # Augment
    if !isfile(assigned_path)
        @info "Augmenting..."
        augmented = augment_table(read_assignments(airr_path), db_dir;
            sequence_type=config.sequence_type)
        write_table_gz(assigned_path, augmented)
    end

    # Filter
    if !isfile(filtered_path)
        filter_table(assigned_path, filtered_path, config.preprocessing_filter;
            stats_path=joinpath(iter_dir, "stats", "filtered.json"))
    end

    # Discover V candidates
    if !isfile(candidates_path)
        table    = read_assignments(filtered_path)
        database = read_fasta_dict(joinpath(db_dir, "V.fasta"))
        write_table(candidates_path, discover_germline(table, database, config))
    end

    # Germline filters
    run_germline_filters(config, iter_dir)

    # J discovery (iteration 1 only)
    if iteration == 1
        j_out = joinpath(iter_dir, "new_J.fasta")
        if !isfile(j_out)
            table     = read_assignments(filtered_path)
            j_database = read_fasta_dict(joinpath(db_dir, "J.fasta"))
            result = discover_j_to_fasta(table, j_database, config, j_out)
            # Fall back to original J if discovery produced nothing
            if isempty(read_fasta(result))
                @info "J discovery produced no candidates; keeping original J database"
                cp(joinpath(db_dir, "J.fasta"), j_out; force=true)
            end
        end
    end
end

function run_germline_filters(config::Config, iter_dir::String)
    candidates_path = joinpath(iter_dir, "candidates.tab")
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
            write_table(joinpath(iter_dir, "new_V_$(prefix)germline.tab"), candidates)
            write_table(joinpath(iter_dir, "annotated_V_$(prefix)germline.tab"), candidates)
            write_fasta(fasta_path, FastaRecord[])
        else
            filtered, annotated = germline_filter!(candidates, criteria; whitelist=wl)
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

    # Choose V source: use discovered germline if non-empty, else fall back to initial DB
    germline_fasta = joinpath(last_iter, "new_V_germline.fasta")
    v_src = if isfile(germline_fasta) && !isempty(read_fasta(germline_fasta))
        germline_fasta
    else
        @info "No new V germline sequences discovered; using original database"
        joinpath("database", "V.fasta")
    end

    cp(v_src, joinpath(final_dir, "database", "V.fasta"); force=true)
    cp(joinpath("database", "D.fasta"), joinpath(final_dir, "database", "D.fasta"); force=true)

    j_src = config.j_discovery.propagate && isfile("iteration-01/new_J.fasta") ?
        "iteration-01/new_J.fasta" : joinpath("database", "J.fasta")
    cp(j_src, joinpath(final_dir, "database", "J.fasta"); force=true)

    db_dir        = joinpath(final_dir, "database")
    airr_path     = joinpath(final_dir, "airr.tsv.gz")
    assigned_path = joinpath(final_dir, "assigned.tsv.gz")
    filtered_path = joinpath(final_dir, "filtered.tsv.gz")

    # Reuse iteration results if we fell back to the original V database
    if v_src == joinpath("database", "V.fasta") && isfile(joinpath(last_iter, "airr.tsv.gz"))
        @info "Reusing iteration results for final (V unchanged)"
        for fname in ("airr.tsv.gz", "assigned.tsv.gz", "filtered.tsv.gz")
            cp(joinpath(last_iter, fname), joinpath(final_dir, fname); force=true)
        end
        src_stats = joinpath(last_iter, "stats", "filtered.json")
        isfile(src_stats) &&
            cp(src_stats, joinpath(final_dir, "stats", "filtered.json"); force=true)
    else
        if !isfile(airr_path)
            run_igblast_on_fasta(db_dir, reads_path, airr_path;
                sequence_type=config.sequence_type, species=config.species,
                penalty=config.mismatch_penalty, threads=Sys.CPU_THREADS,
                limit=config.limit)
        end
        if !isfile(assigned_path)
            write_table_gz(assigned_path,
                augment_table(read_assignments(airr_path), db_dir;
                    sequence_type=config.sequence_type))
        end
        if !isfile(filtered_path)
            filter_table(assigned_path, filtered_path, config.preprocessing_filter;
                stats_path=joinpath(final_dir, "stats", "filtered.json"))
        end
    end

    @info "Final output in $final_dir/"
end
