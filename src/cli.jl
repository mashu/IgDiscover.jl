# Command-line interface and PackageCompiler app entry point.

function argparse_settings()
    s = ArgParse.ArgParseSettings(
        "IgDiscover.jl — antibody repertoire analysis and V gene discovery",
        version = "IgDiscover.jl 0.3.0",
        add_version = false,
    )
    s.commands_are_required = false
    s.epilog = """
        Configuration:
          Edit igdiscover.toml in the analysis directory before running.
          Set limit=1000 for a quick test run, limit=0 for the full dataset.
        """
    @ArgParse.add_arg_table! s begin
        "init"
            action = :command
            help = "Initialize an analysis directory"
        "run"
            action = :command
            help = "Run the pipeline"
        "clonotypes"
            action = :command
            help = "Compute clonotypes from a filtered table"
        "clonoquery"
            action = :command
            help = "Query a reference table by clonotype"
        "count"
            action = :command
            help = "Compute gene expression counts"
        "dbdiff"
            action = :command
            help = "Compare two FASTA databases"
        "upstream"
            action = :command
            help = "Compute UTR/leader consensus per V gene"
        "plotalleles"
            action = :command
            help = "Compute allele co-expression matrix"
        "haplotype"
            action = :command
            help = "Determine haplotypes from allele co-occurrences"
        "clusterplot"
            action = :command
            help = "Cluster sequences per gene and show distance heatmap"
        "dereplicate"
            action = :command
            help = "Deduplicate sequences and remove barcodes"
        "merge"
            action = :command
            help = "Merge paired-end reads (PEAR wrapper)"
        "union"
            action = :command
            help = "Compute union of sequences in multiple FASTA files"
        "version"
            action = :command
            help = "Show version"
    end

    @ArgParse.add_arg_table! s["init"] begin
        "dir"
            help = "Analysis directory to create"
            required = true
        "database"
            help = "Path to database directory"
            required = true
        "reads"
            help = "Path to reads file"
            required = true
    end

    @ArgParse.add_arg_table! s["run"] begin
        "dir"
            help = "Analysis directory (default: .)"
            required = false
            default = "."
    end

    @ArgParse.add_arg_table! s["clonotypes"] begin
        "--sort"
            help = "Sort by group size (largest first)"
            action = :store_true
        "--limit"
            help = "Print only the first N clonotypes"
            arg_type = Int
            default = 0
        "--mismatches"
            help = "Max CDR3 mismatches"
            arg_type = Float64
            default = 1.0
        "--v-shm-threshold"
            help = "V SHM threshold for _mindiffrate"
            arg_type = Float64
            default = 5.0
        "--members"
            help = "Write member table to FILE"
            arg_type = String
            default = ""
        "table"
            help = "Input filtered table"
            required = true
        "output"
            help = "Output clonotype table"
            required = true
    end

    @ArgParse.add_arg_table! s["clonoquery"] begin
        "--minimum-count"
            help = "Discard ref rows with count < N"
            arg_type = Int
            default = 1
        "--mismatches"
            help = "Max CDR3 mismatches"
            arg_type = Float64
            default = 1.0
        "--aa"
            help = "Compare CDR3 at amino-acid level"
            action = :store_true
        "--summary"
            help = "Write summary to FILE"
            arg_type = String
            default = ""
        "reftable"
            help = "Reference table"
            required = true
        "querytable"
            help = "Query table"
            required = true
    end

    @ArgParse.add_arg_table! s["count"] begin
        "--gene"
            help = "Gene type: V, D, or J"
            arg_type = String
            default = "V"
        "--database"
            help = "FASTA file with gene names"
            arg_type = String
            default = ""
        "--d-evalue"
            help = "Max D gene E-value"
            arg_type = Float64
            default = NaN
        "--d-coverage"
            help = "Min D coverage"
            arg_type = Float64
            default = NaN
        "--allele-ratio"
            help = "Min allele expression ratio"
            arg_type = Float64
            default = NaN
        "table"
            help = "Input table"
            required = true
    end

    @ArgParse.add_arg_table! s["dbdiff"] begin
        "a"
            help = "First FASTA file"
            required = true
        "b"
            help = "Second FASTA file"
            required = true
    end

    @ArgParse.add_arg_table! s["upstream"] begin
        "--max-v-errors"
            help = "Allow PERCENT errors in V"
            arg_type = Float64
            default = 1.0
        "--min-consensus-size"
            help = "Min sequences for consensus"
            arg_type = Int
            default = 1
        "--consensus-threshold"
            help = "Consensus threshold (%%)"
            arg_type = Float64
            default = 75.0
        "--no-ambiguous"
            help = "Discard consensus with N bases"
            action = :store_true
        "--part"
            help = "UTR, leader, or UTR+leader"
            arg_type = String
            default = "UTR+leader"
        "table"
            help = "Input table"
            required = true
        "output"
            help = "Output FASTA"
            required = true
    end

    @ArgParse.add_arg_table! s["plotalleles"] begin
        "--x"
            help = "Gene type on x axis"
            arg_type = String
            default = "V"
        "--gene"
            help = "Gene type on y axis"
            arg_type = String
            default = "J"
        "--d-evalue"
            help = "Max D E-value"
            arg_type = Float64
            default = 1e-4
        "--d-coverage"
            help = "Min D coverage"
            arg_type = Float64
            default = 65.0
        "--database"
            help = "Restrict x-axis to this FASTA"
            arg_type = String
            default = ""
        "--order"
            help = "Sort x-axis by this FASTA"
            arg_type = String
            default = ""
        "alleles"
            help = "Comma-separated alleles for y axis"
            required = true
        "table"
            help = "Input table"
            required = true
        "output"
            help = "Output TSV"
            required = true
    end

    @ArgParse.add_arg_table! s["haplotype"] begin
        "--v-gene"
            help = "V gene to use for haplotyping"
            arg_type = String
            default = ""
        "--d-evalue"
            help = "Max D E-value"
            arg_type = Float64
            default = 1e-4
        "--d-coverage"
            help = "Min D coverage"
            arg_type = Float64
            default = 65.0
        "--restrict"
            help = "Restrict to genes in this FASTA"
            arg_type = String
            default = ""
        "--order"
            help = "Sort by record order in FASTA"
            arg_type = String
            default = ""
        "table"
            help = "Input filtered table"
            required = true
        "output"
            help = "Output TSV"
            required = true
    end

    @ArgParse.add_arg_table! s["clusterplot"] begin
        "--minimum-group-size"
            help = "Skip genes with fewer sequences"
            arg_type = Int
            default = 200
        "--size"
            help = "Max sequences to display"
            arg_type = Int
            default = 300
        "--type"
            help = "Gene type: V, D, or J"
            arg_type = String
            default = "V"
        "--gene"
            help = "Plot only this gene (repeatable via comma)"
            arg_type = String
            default = ""
        "--ignore-j"
            help = "Include rows without J assignment"
            action = :store_true
        "table"
            help = "Input table"
            required = true
    end

    @ArgParse.add_arg_table! s["dereplicate"] begin
        "--limit"
            help = "Process only first N reads"
            arg_type = Int
            default = 0
        "--trim-g"
            help = "Trim leading G nucleotides (RACE artifact)"
            action = :store_true
        "--minimum-length"
            help = "Minimum sequence length"
            arg_type = Int
            default = 0
        "--barcode-length"
            help = "Barcode length (positive=5', negative=3')"
            arg_type = Int
            default = 0
        "--json"
            help = "Write statistics to JSON file"
            arg_type = String
            default = ""
        "fastx"
            help = "Input FASTA/FASTQ file"
            required = true
    end

    @ArgParse.add_arg_table! s["merge"] begin
        "--threads"
            help = "Number of threads"
            arg_type = Int
            default = 0
        "reads1"
            help = "Forward reads FASTQ file"
            required = true
        "reads2"
            help = "Reverse reads FASTQ file"
            required = true
        "output"
            help = "Output merged FASTQ.gz"
            required = true
    end

    @ArgParse.add_arg_table! s["union"] begin
        "fasta"
            help = "FASTA files to merge"
            nargs = '+'
            required = true
    end

    return s
end

# ─── Command handlers (dispatch by name via Val) ───

handle_command(::Val{:init}, parsed) = begin
    a = parsed["init"]
    isdir(a["database"]) || error("Database directory not found: $(a["database"])")
    isfile(a["reads"]) || error("Reads file not found: $(a["reads"])")
    init_analysis(a["dir"], a["database"], a["reads"])
end

handle_command(::Val{:run}, parsed) = begin
    dir = parsed["run"]["dir"]
    isdir(dir) || error("Analysis directory not found: $dir")
    run_pipeline(dir)
end

handle_command(::Val{:version}, _) = println("IgDiscover.jl 0.3.0")

handle_command(::Val{:clonotypes}, parsed) = begin
    a = parsed["clonotypes"]
    table = read_assignments(a["table"])
    caller = ClonotypeCaller(
        max_mismatches=round(Int, a["mismatches"]),
        v_shm_threshold=a["v-shm-threshold"],
        sort_by_size=a["sort"],
    )
    ct_df, clusters = call_clonotypes(table, caller)
    limit = a["limit"]
    if limit > 0
        ct_df = ct_df[1:min(limit, nrow(ct_df)), :]
        clusters = clusters[1:min(limit, length(clusters))]
    end
    write_table(a["output"], ct_df)
    @info "Wrote $(nrow(ct_df)) clonotypes to $(a["output"])"
    members_path = a["members"]
    !isempty(members_path) && (write_members(members_path, table, clusters);
        @info "Wrote members to $members_path")
end

handle_command(::Val{:clonoquery}, parsed) = begin
    a = parsed["clonoquery"]
    cq = ClonoQuery(mismatches=a["mismatches"], use_aa=a["aa"], min_count=a["minimum-count"])
    results = cq(read_assignments(a["querytable"]), read_assignments(a["reftable"]))
    write_clonoquery(stdout, results)
    summary_path = a["summary"]
    if !isempty(summary_path)
        write_table(summary_path, clonoquery_summary(results, read_assignments(a["reftable"])))
        @info "Wrote summary to $summary_path"
    end
end

handle_command(::Val{:count}, parsed) = begin
    a = parsed["count"]
    gene_sym = Symbol(uppercase(a["gene"]))
    gene_names = !isempty(a["database"]) ?
        sort([r.name for r in read_fasta(a["database"])]; by=natural_sort_key) : nothing
    counter = ExpressionCounter(gene=gene_sym,
        d_evalue=isnan(a["d-evalue"]) ? nothing : a["d-evalue"],
        d_coverage=isnan(a["d-coverage"]) ? nothing : a["d-coverage"],
        allele_ratio=isnan(a["allele-ratio"]) ? nothing : a["allele-ratio"])
    counts = counter(read_assignments(a["table"]); gene_names=gene_names)
    for i in 1:nrow(counts)
        println(counts.gene[i], '\t', counts.count[i])
    end
end

handle_command(::Val{:dbdiff}, parsed) = begin
    a = parsed["dbdiff"]
    a_recs = read_fasta(a["a"])
    b_recs = read_fasta(a["b"])
    for (recs, path) in ((a_recs, a["a"]), (b_recs, a["b"]))
        for n in check_duplicate_names(recs); println("Duplicate name in $path: $n"); end
        for (n, o) in check_duplicate_sequences(recs); println("Duplicate seq in $path: $n = $o"); end
    end
    diff = DatabaseComparator()(a_recs, b_recs)
    format_diff(stdout, diff)
end

handle_command(::Val{:upstream}, parsed) = begin
    a = parsed["upstream"]
    part = a["part"] == "UTR" ? :UTR : a["part"] == "leader" ? :leader : :UTR_leader
    analyzer = UpstreamAnalyzer(UpstreamParams(
        max_v_errors=a["max-v-errors"], min_consensus_size=a["min-consensus-size"],
        consensus_threshold=a["consensus-threshold"], discard_ambiguous=a["no-ambiguous"],
        part=part))
    records = analyzer(read_assignments(a["table"]))
    write_fasta(a["output"], records)
    @info "Wrote $(length(records)) consensus sequences to $(a["output"])"
end

handle_command(::Val{:plotalleles}, parsed) = begin
    a = parsed["plotalleles"]
    au = AlleleUsage(x_gene=Symbol(uppercase(a["x"])), y_gene=Symbol(uppercase(a["gene"])),
        filter_criteria=AlleleFilterCriteria(d_evalue=a["d-evalue"], d_coverage=a["d-coverage"]))
    db_names = !isempty(a["database"]) ? [r.name for r in read_fasta(a["database"])] : nothing
    order_names = !isempty(a["order"]) ? [first(split(r.name, '*')) for r in read_fasta(a["order"])] : nothing
    matrix = au(read_assignments(a["table"]), String.(split(a["alleles"], ','));
                database_names=db_names, order_names=order_names)
    write_table(a["output"], matrix)
    @info "Wrote allele matrix to $(a["output"])"
end

handle_command(::Val{:haplotype}, parsed) = begin
    a = parsed["haplotype"]
    ha = HaplotypeAnalyzer(d_evalue=a["d-evalue"], d_coverage=a["d-coverage"])
    restrict = !isempty(a["restrict"]) ? Set(r.name for r in read_fasta(a["restrict"])) : nothing
    order = !isempty(a["order"]) ? [r.name for r in read_fasta(a["order"])] : nothing
    v_gene = !isempty(a["v-gene"]) ? a["v-gene"] : nothing
    blocks = ha(read_assignments(a["table"]); restrict_names=restrict,
                gene_order=order, v_gene=v_gene)
    open(a["output"], "w") do io
        for (i, block) in enumerate(blocks)
            print(io, format_tsv(block; header=(i == 1)))
        end
    end
    @info "Wrote haplotype table to $(a["output"])"
end

handle_command(::Val{:clusterplot}, parsed) = begin
    a = parsed["clusterplot"]
    gene_type = Symbol(uppercase(a["type"]))
    genes = isempty(a["gene"]) ? nothing : String.(split(a["gene"], ','))
    plotter = ClusterPlotter(min_group_size=a["minimum-group-size"],
        max_display=a["size"], gene_type=gene_type, ignore_j=a["ignore-j"])
    results = plotter(read_assignments(a["table"]); genes=genes)
    for r in results
        println(render_heatmap(r))
    end
    @info "Processed $(length(results)) genes"
end

handle_command(::Val{:dereplicate}, parsed) = begin
    a = parsed["dereplicate"]
    derep = Dereplicator(DereplicateParams(
        barcode_length=a["barcode-length"], trim_g=a["trim-g"],
        minimum_length=a["minimum-length"], limit=a["limit"]))
    records, stats = derep(a["fastx"])
    write_dereplicated(stdout, records, a["barcode-length"] != 0)
    @info "$(stats.total_reads) reads → $(stats.unique_sequences) unique ($(stats.too_short) too short)"
    json_path = a["json"]
    if !isempty(json_path)
        open(json_path, "w") do io
            JSON3.pretty(io, Dict("groups_written" => stats.unique_sequences))
            println(io)
        end
    end
end

handle_command(::Val{:merge}, parsed) = begin
    a = parsed["merge"]
    merger = ReadMerger(program=:pear, threads=a["threads"])
    stats = merger(a["reads1"], a["reads2"], a["output"])
    @info "Merged $(stats.merged) of $(stats.total) read pairs"
end

handle_command(::Val{:union}, parsed) = begin
    a = parsed["union"]
    merged = SequenceUnion()(a["fasta"])
    for rec in merged
        println(">$(rec.name)\n$(rec.sequence)")
    end
end

# ─── Entry point ───

function cli_main()
    s = argparse_settings()
    parsed = ArgParse.parse_args(ARGS, s)
    cmd = parsed["%COMMAND%"]

    if cmd === nothing
        ArgParse.show_help(s, exit_when_done = false)
        return 0
    end

    handle_command(Val(Symbol(cmd)), parsed)
    return 0
end

"""
    julia_main()::Cint

Entry point for PackageCompiler standalone executable.
"""
function julia_main()::Cint
    try
        return Cint(cli_main())
    catch e
        Base.invokelatest(Base.display_error, e, catch_backtrace())
        return Cint(1)
    end
end
