#!/usr/bin/env julia
# IgDiscover.jl command-line interface
#
# Usage:
#   julia --project=. bin/igdiscover.jl init <dir> <database> <reads>
#   julia --project=. bin/igdiscover.jl run [<dir>]

using Pkg
script_dir = @__DIR__
project_dir = dirname(script_dir)
Pkg.activate(project_dir)

using IgDiscover

function print_usage()
    println("""
    IgDiscover.jl — antibody repertoire analysis and V gene discovery

    Commands:
      init <dir> <database> <reads>   Initialize an analysis directory
      run  [<dir>]                    Run the pipeline (default: current directory)
      version                         Show version

    Examples:
      julia bin/igdiscover.jl init my_analysis database/ reads.fasta.gz
      julia bin/igdiscover.jl run my_analysis

    Configuration:
      Edit igdiscover.toml in the analysis directory before running.
      Set limit=1000 for a quick test run, limit=0 for the full dataset.
    """)
end

function main()
    isempty(ARGS) && (print_usage(); exit(0))

    cmd = ARGS[1]

    if cmd == "init"
        length(ARGS) < 4 && (println("Usage: igdiscover init <dir> <database_dir> <reads>"); exit(1))
        dir, db, reads = ARGS[2], ARGS[3], ARGS[4]
        isdir(db) || error("Database directory not found: $db")
        isfile(reads) || error("Reads file not found: $reads")
        init_analysis(dir, db, reads)

    elseif cmd == "run"
        dir = length(ARGS) >= 2 ? ARGS[2] : "."
        isdir(dir) || error("Analysis directory not found: $dir")
        run_pipeline(dir)

    elseif cmd == "version"
        println("IgDiscover.jl 0.1.0")

    elseif cmd in ("-h", "--help", "help")
        print_usage()

    else
        println("Unknown command: $cmd")
        print_usage()
        exit(1)
    end
end

main()
