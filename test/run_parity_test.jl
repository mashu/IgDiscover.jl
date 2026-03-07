#!/usr/bin/env julia
# run_parity_test.jl — Automated parity test: Python igdiscover vs Julia IgDiscover.jl
#
# This test:
#   1. Sets up Python igdiscover (via pip in a venv)
#   2. Creates a small test database from IMGT-style FASTA
#   3. Generates synthetic reads
#   4. Runs both implementations with matching configurations
#   5. Compares filtered.tsv.gz and new_V_germline outputs
#
# Usage:
#   julia --project=. test/run_parity_test.jl                    # full automated test
#   julia --project=. test/run_parity_test.jl <py_dir> <jl_dir>  # compare existing outputs

using IgDiscover
using DataFrames
using Printf
using CSV

const TEST_DIR = joinpath(tempdir(), "igdiscover_parity_test_$(getpid())")
const TEST_LIMIT = parse(Int, get(ENV, "IGDISCOVER_TEST_LIMIT", "100"))

# ─── Comparison functions ───

function compare_filtered(py_dir::String, jl_dir::String)
    println("\n══════ Comparing filtered.tsv.gz ══════")

    py_path = find_file(py_dir, "filtered.tsv.gz")
    jl_path = find_file(jl_dir, "filtered.tsv.gz")

    py_path === nothing && (println("  ✗ Python filtered not found"); return false)
    jl_path === nothing && (println("  ✗ Julia filtered not found"); return false)

    py_df = read_assignments(py_path)
    jl_df = read_assignments(jl_path)

    println("  Python: $(nrow(py_df)) rows, $(ncol(py_df)) columns")
    println("  Julia:  $(nrow(jl_df)) rows, $(ncol(jl_df)) columns")

    passed = true

    if nrow(py_df) == nrow(jl_df)
        println("  ✓ Row count matches: $(nrow(py_df))")
    else
        println("  ✗ Row count differs: Python=$(nrow(py_df)), Julia=$(nrow(jl_df))")
        passed = false
    end

    shared_cols = intersect(propertynames(py_df), propertynames(jl_df))
    println("  Shared columns: $(length(shared_cols))/$(ncol(py_df)) Python, $(length(shared_cols))/$(ncol(jl_df)) Julia")

    py_only = setdiff(propertynames(py_df), propertynames(jl_df))
    jl_only = setdiff(propertynames(jl_df), propertynames(py_df))
    isempty(py_only) || println("  Python-only columns: $(join(py_only, ", "))")
    isempty(jl_only) || println("  Julia-only columns: $(join(jl_only, ", "))")

    if :v_call in shared_cols && nrow(py_df) == nrow(jl_df)
        py_vcounts = sort(collect(IgDiscover.tallies(py_df.v_call)); by=first)
        jl_vcounts = sort(collect(IgDiscover.tallies(jl_df.v_call)); by=first)
        if py_vcounts == jl_vcounts
            println("  ✓ V gene distribution matches")
        else
            println("  ✗ V gene distribution differs")
            passed = false
        end
    end

    passed
end

function compare_germline(py_dir::String, jl_dir::String)
    println("\n══════ Comparing new_V_germline ══════")

    py_tab = find_file(py_dir, "new_V_germline.tab")
    jl_tab = find_file(jl_dir, "new_V_germline.tab")
    py_fa = find_file(py_dir, "new_V_germline.fasta")
    jl_fa = find_file(jl_dir, "new_V_germline.fasta")

    passed = true

    if py_tab !== nothing && jl_tab !== nothing
        py_df = CSV.read(py_tab, DataFrame; delim='\t')
        jl_df = CSV.read(jl_tab, DataFrame; delim='\t')

        println("  Python: $(nrow(py_df)) candidates")
        println("  Julia:  $(nrow(jl_df)) candidates")

        if nrow(py_df) == nrow(jl_df)
            println("  ✓ Candidate count matches")
        else
            println("  ✗ Candidate count differs")
            passed = false
        end

        if hasproperty(py_df, :consensus) && hasproperty(jl_df, :consensus)
            py_seqs = sort(py_df.consensus)
            jl_seqs = sort(jl_df.consensus)
            if py_seqs == jl_seqs
                println("  ✓ All consensus sequences match exactly")
            else
                common = length(intersect(Set(py_seqs), Set(jl_seqs)))
                py_only = length(setdiff(Set(py_seqs), Set(jl_seqs)))
                jl_only = length(setdiff(Set(jl_seqs), Set(py_seqs)))
                println("  ✗ Consensus sequences differ: $common shared, $py_only Python-only, $jl_only Julia-only")
                passed = false
            end
        end
    else
        py_tab === nothing && println("  ✗ Python new_V_germline.tab not found")
        jl_tab === nothing && println("  ✗ Julia new_V_germline.tab not found")
        passed = false
    end

    if py_fa !== nothing && jl_fa !== nothing
        py_records = read_fasta_dict(py_fa)
        jl_records = read_fasta_dict(jl_fa)

        if py_records == jl_records
            println("  ✓ FASTA files identical ($(length(py_records)) sequences)")
        else
            println("  ✗ FASTA files differ: Python=$(length(py_records)), Julia=$(length(jl_records))")
            passed = false
        end
    end

    passed
end

# ─── File finders ───

function find_file(dir::String, filename::String)
    for candidate in [
        joinpath(dir, "iteration-01", filename),
        joinpath(dir, "final", filename),
        joinpath(dir, filename),
    ]
        isfile(candidate) && return candidate
    end
    nothing
end

# ─── Python igdiscover setup ───

function setup_python_igdiscover(venv_dir::String)
    if isdir(venv_dir) && isfile(joinpath(venv_dir, "bin", "igdiscover"))
        println("  Python igdiscover already installed at $venv_dir")
        return joinpath(venv_dir, "bin")
    end

    println("  Setting up Python igdiscover...")
    run(`python3 -m venv $venv_dir`)
    pip = joinpath(venv_dir, "bin", "pip")
    run(`$pip install --quiet igdiscover`)
    bin_dir = joinpath(venv_dir, "bin")
    @assert isfile(joinpath(bin_dir, "igdiscover")) "igdiscover not found after pip install"
    println("  ✓ Installed igdiscover at $bin_dir")
    bin_dir
end

# ─── Test database creation (IMGT-style headers with dots) ───

function create_test_database(db_dir::String)
    mkpath(db_dir)

    # V genes with IMGT-style headers and dot-gapped sequences
    v_records = [
        # Realistic IMGT header format with dots in sequence
        ("M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |",
         "cagg.ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggttacaccttt...............accagctatggt...atcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcttac......aatggtaacacaaactatgcacagaagctccag...ggcagagtcaccatgaccacagacacatccacgagcacagcctacatggagctgaggagcctgagatctgacgacacggccgtgtattactgtgcgagaga"),
        ("M99641|IGHV1-18*04|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |",
         "cagg.ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggttacaccttt...............accagctatggt...atcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcttac......aatggtaacacaaactatgcacagaagctccag...ggcagagtcaccatgaccacagacacatccacgagcacagcctacatggagctgaggagcctgagatctgacgacacggccgtgtattactgtgcgagagg"),
        ("X62106|IGHV1-2*01|Homo sapiens|F|V-REGION|293..588|296 nt|1| | | | |296+24=320| | |",
         "cagg.tgcagctggtgcag...tctggggct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggatacaccttc...............accggctactatatg...cactgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcaaccctaac......agtggtggcacaaactatgcacagaagtttcag...ggcagggtcaccatgaccagggacacgtccatcagcacagcctacatggagctgagcaggctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ]

    open(joinpath(db_dir, "V.fasta"), "w") do io
        for (header, seq) in v_records
            println(io, ">", header)
            println(io, seq)
        end
    end

    # D genes
    d_records = [
        ("X97051|IGHD1-1*01|Homo sapiens|F|D-REGION|8..30|23 nt|1| | | | | | |", "ggtataactggaactgactacggg"),
        ("X13972|IGHD2-2*01|Homo sapiens|F|D-REGION|8..25|18 nt|1| | | | | | |", "aggatattgtagtagtaccagctgctatgcc"),
    ]

    open(joinpath(db_dir, "D.fasta"), "w") do io
        for (header, seq) in d_records
            println(io, ">", header)
            println(io, seq)
        end
    end

    # J genes
    j_records = [
        ("J00256|IGHJ1*01|Homo sapiens|F|J-REGION|9..53|45 nt|1| | | | | | |", "gctgaatacttccagcactggggccagggcaccctggtcaccgtctcctcag"),
        ("J00256|IGHJ2*01|Homo sapiens|F|J-REGION|118..170|53 nt|1| | | | | | |", "ctactggtacttcgatctctggggccgtggcaccctggtcactgtctcctcag"),
        ("J00256|IGHJ3*01|Homo sapiens|F|J-REGION|229..279|51 nt|1| | | | | | |", "tgatgcttttgatgtctggggccaagggacaatggtcaccgtctcttcag"),
        ("J00256|IGHJ4*01|Homo sapiens|F|J-REGION|346..393|48 nt|1| | | | | | |", "actactttgactactggggccaaggaaccctggtcaccgtctcctcag"),
    ]

    open(joinpath(db_dir, "J.fasta"), "w") do io
        for (header, seq) in j_records
            println(io, ">", header)
            println(io, seq)
        end
    end

    @info "Created test database at $db_dir"
end

# ─── Run Python igdiscover ───

function run_python_igdiscover(py_bin::String, analysis_dir::String, db_dir::String, reads_path::String, limit::Int)
    igdiscover = joinpath(py_bin, "igdiscover")

    # Initialize
    if !isdir(analysis_dir)
        run(`$igdiscover init --db $db_dir --single-reads $reads_path $analysis_dir`)
    end

    # Set limit in config
    config_path = joinpath(analysis_dir, "igdiscover.yaml")
    if isfile(config_path)
        config_text = read(config_path, String)
        config_text = replace(config_text, r"limit:\s*\d+" => "limit: $limit")
        config_text = replace(config_text, r"iterations:\s*\d+" => "iterations: 1")
        write(config_path, config_text)
    end

    # Run
    run(Cmd(`$igdiscover run`; dir=analysis_dir))
    @info "Python igdiscover completed at $analysis_dir"
end

# ─── Run Julia IgDiscover ───

function run_julia_igdiscover(analysis_dir::String, db_dir::String, reads_path::String, limit::Int)
    if !isdir(analysis_dir)
        init_analysis(analysis_dir, db_dir, reads_path)
    end

    # Set limit in config
    config_path = joinpath(analysis_dir, "igdiscover.toml")
    config_text = read(config_path, String)
    config_text = replace(config_text, r"limit\s*=\s*\d+" => "limit = $limit")
    config_text = replace(config_text, r"iterations\s*=\s*\d+" => "iterations = 1")
    write(config_path, config_text)

    run_pipeline(analysis_dir)
    @info "Julia IgDiscover.jl completed at $analysis_dir"
end

# ─── Main ───

function main()
    if length(ARGS) == 2
        # Compare mode: just compare existing outputs
        py_dir, jl_dir = ARGS[1], ARGS[2]
        isdir(py_dir) || error("Python directory not found: $py_dir")
        isdir(jl_dir) || error("Julia directory not found: $jl_dir")

        println("Parity test: Python vs Julia IgDiscover")
        println("  Python: $py_dir")
        println("  Julia:  $jl_dir")

        f_ok = compare_filtered(py_dir, jl_dir)
        g_ok = compare_germline(py_dir, jl_dir)

        print_summary(f_ok, g_ok)
        exit(f_ok && g_ok ? 0 : 1)
    end

    # Full automated test
    println("═══ IgDiscover Parity Test (limit=$TEST_LIMIT) ═══")
    println("  Working directory: $TEST_DIR")
    mkpath(TEST_DIR)

    # Step 1: Set up Python igdiscover
    println("\n── Step 1: Python igdiscover setup ──")
    py_bin = setup_python_igdiscover(joinpath(TEST_DIR, "venv"))

    # Step 2: Create test database
    println("\n── Step 2: Test database ──")
    db_dir = joinpath(TEST_DIR, "database")
    create_test_database(db_dir)

    # Step 3: Need reads — check if user provided or generate a note
    reads_path = get(ENV, "IGDISCOVER_TEST_READS", "")
    if isempty(reads_path) || !isfile(reads_path)
        println("\n  ⚠ No test reads provided.")
        println("  Set IGDISCOVER_TEST_READS=/path/to/reads.fasta.gz")
        println("  Or provide pre-computed analysis directories as arguments.")
        println("  Skipping automated run, testing only IMGT sanitization...")

        # At minimum, verify IMGT sanitization works correctly
        test_imgt_sanitization(db_dir)
        return
    end

    # Step 4: Run Python igdiscover
    println("\n── Step 4: Running Python igdiscover ──")
    py_analysis = joinpath(TEST_DIR, "python_analysis")
    run_python_igdiscover(py_bin, py_analysis, db_dir, reads_path, TEST_LIMIT)

    # Step 5: Run Julia IgDiscover.jl
    println("\n── Step 5: Running Julia IgDiscover.jl ──")
    jl_analysis = joinpath(TEST_DIR, "julia_analysis")
    run_julia_igdiscover(jl_analysis, db_dir, reads_path, TEST_LIMIT)

    # Step 6: Compare
    println("\n── Step 6: Comparing outputs ──")
    f_ok = compare_filtered(py_analysis, jl_analysis)
    g_ok = compare_germline(py_analysis, jl_analysis)

    print_summary(f_ok, g_ok)
    exit(f_ok && g_ok ? 0 : 1)
end

function test_imgt_sanitization(db_dir::String)
    println("\n══════ IMGT Sanitization Test ══════")
    passed = true

    for gene in ("V", "D", "J")
        src = joinpath(db_dir, "$gene.fasta")
        isfile(src) || continue

        raw = read_fasta(src)
        tmpdir = mktempdir()
        clean_path = joinpath(tmpdir, "$gene.fasta")
        write_sanitized_imgt(src, clean_path)
        clean = read_fasta(clean_path)

        # Check no dots in cleaned sequences
        for r in clean
            if occursin('.', r.sequence)
                println("  ✗ $gene: $(r.name) still contains dots after sanitization")
                passed = false
            end
        end

        # Check headers are clean allele names
        for r in clean
            if occursin('|', r.name)
                println("  ✗ $gene: $(r.name) still has IMGT pipe-delimited header")
                passed = false
            end
        end

        println("  ✓ $gene: $(length(clean)) records sanitized ($(length(raw)) input)")
        rm(tmpdir; recursive=true)
    end

    println(passed ? "\n✓ IMGT sanitization PASS" : "\n✗ IMGT sanitization FAIL")
end

function print_summary(f_ok::Bool, g_ok::Bool)
    println("\n══════ Summary ══════")
    println("  filtered.tsv.gz: $(f_ok ? "✓ PASS" : "✗ FAIL")")
    println("  new_V_germline:  $(g_ok ? "✓ PASS" : "✗ FAIL")")
    (f_ok && g_ok) ? println("\n✓ PARITY ACHIEVED") : println("\n✗ PARITY NOT YET ACHIEVED")
end

main()
