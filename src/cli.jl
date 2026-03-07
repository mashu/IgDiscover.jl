# Command-line interface and PackageCompiler app entry point.
# The compiled standalone app calls IgDiscover.julia_main()::Cint.

function _argparse_settings()
    s = ArgParse.ArgParseSettings(
        "IgDiscover.jl — antibody repertoire analysis and V gene discovery",
        version = "IgDiscover.jl 0.1.0",
        add_version = false,  # we use the "version" command instead
    )
    s.commands_are_required = false  # no args → show usage instead of error
    s.epilog = """
        Configuration:
          Edit igdiscover.toml in the analysis directory before running.
          Set limit=1000 for a quick test run, limit=0 for the full dataset.
        """
    @ArgParse.add_arg_table! s begin
        "init"
            action = :command
            help = "Initialize an analysis directory with database, reads, and default config"
        "run"
            action = :command
            help = "Run the pipeline (default: current directory)"
        "version"
            action = :command
            help = "Show version"
    end

    # init: dir, database_dir, reads
    @ArgParse.add_arg_table! s["init"] begin
        "dir"
            help = "Analysis directory to create"
            required = true
        "database"
            help = "Path to database directory (V.fasta, D.fasta, J.fasta)"
            required = true
        "reads"
            help = "Path to reads file (FASTA/FASTQ, optionally gzipped)"
            required = true
    end

    # run: optional dir (default ".")
    @ArgParse.add_arg_table! s["run"] begin
        "dir"
            help = "Analysis directory (default: current directory)"
            required = false
            default = "."
    end

    return s
end

function _run_init(args)
    dir = args["init"]["dir"]
    db = args["init"]["database"]
    reads = args["init"]["reads"]
    isdir(db) || error("Database directory not found: $db")
    isfile(reads) || error("Reads file not found: $reads")
    init_analysis(dir, db, reads)
end

function _run_run(args)
    dir = args["run"]["dir"]
    isdir(dir) || error("Analysis directory not found: $dir")
    run_pipeline(dir)
end

function _run_version()
    println("IgDiscover.jl 0.1.0")
end

function _cli_main()
    s = _argparse_settings()
    parsed = ArgParse.parse_args(ARGS, s)
    cmd = parsed["%COMMAND%"]

    if cmd === nothing
        ArgParse.show_help(s, exit_when_done = false)
        return 0
    end

    if cmd == "init"
        _run_init(parsed)
    elseif cmd == "run"
        _run_run(parsed)
    elseif cmd == "version"
        _run_version()
    else
        # Should not happen if ArgParse is configured correctly
        println("Unknown command: $cmd")
        return 1
    end
    return 0
end

"""
    julia_main()::Cint

Entry point for the PackageCompiler-built standalone executable.
Parses ARGS (set by the C driver), dispatches to init/run/version, and returns an exit code.
"""
function julia_main()::Cint
    try
        return Cint(_cli_main())
    catch e
        Base.invokelatest(Base.display_error, e, catch_backtrace())
        return Cint(1)
    end
end
