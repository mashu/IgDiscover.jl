#!/usr/bin/env julia
# run_parity_test.jl — Compare Python igdiscover and Julia IgDiscover.jl outputs
#
# Usage:
#   julia --project=. test/run_parity_test.jl /path/to/python/analysis_dir /path/to/julia/analysis_dir

using IgDiscover
using DataFrames
using Printf
using CSV

function compare_filtered(py_dir, jl_dir)
    println("\n══════ Comparing filtered.tsv.gz ══════")

    py_path = find_filtered(py_dir)
    jl_path = find_filtered(jl_dir)

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
        py_vcounts = sort(collect(IgDiscover.tallies(py_df.v_call)); by = first)
        jl_vcounts = sort(collect(IgDiscover.tallies(jl_df.v_call)); by = first)
        if py_vcounts == jl_vcounts
            println("  ✓ V gene distribution matches")
        else
            println("  ✗ V gene distribution differs")
            passed = false
        end
    end

    passed
end

function compare_germline(py_dir, jl_dir)
    println("\n══════ Comparing new_V_germline ══════")

    py_tab = find_germline_tab(py_dir)
    jl_tab = find_germline_tab(jl_dir)
    py_fa = find_germline_fasta(py_dir)
    jl_fa = find_germline_fasta(jl_dir)

    passed = true

    if py_tab !== nothing && jl_tab !== nothing
        py_df = CSV.read(py_tab, DataFrame; delim = '\t')
        jl_df = CSV.read(jl_tab, DataFrame; delim = '\t')

        println("  Python: $(nrow(py_df)) candidates")
        println("  Julia:  $(nrow(jl_df)) candidates")

        if nrow(py_df) == nrow(jl_df)
            println("  ✓ Candidate count matches")
        else
            println("  ✗ Candidate count differs")
            passed = false
        end

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

function find_filtered(dir)
    for candidate in [
        joinpath(dir, "final", "filtered.tsv.gz"),
        joinpath(dir, "iteration-01", "filtered.tsv.gz"),
    ]
        isfile(candidate) && return candidate
    end
    nothing
end

function find_germline_tab(dir)
    for candidate in [
        joinpath(dir, "iteration-01", "new_V_germline.tab"),
        joinpath(dir, "final", "new_V_germline.tab"),
    ]
        isfile(candidate) && return candidate
    end
    nothing
end

function find_germline_fasta(dir)
    for candidate in [
        joinpath(dir, "iteration-01", "new_V_germline.fasta"),
        joinpath(dir, "final", "database", "V.fasta"),
    ]
        isfile(candidate) && return candidate
    end
    nothing
end

function main()
    if length(ARGS) != 2
        println("Usage: julia --project=. test/run_parity_test.jl <python_dir> <julia_dir>")
        exit(1)
    end

    py_dir, jl_dir = ARGS[1], ARGS[2]
    isdir(py_dir) || error("Python directory not found: $py_dir")
    isdir(jl_dir) || error("Julia directory not found: $jl_dir")

    println("Parity test: Python vs Julia IgDiscover")
    println("  Python: $py_dir")
    println("  Julia:  $jl_dir")

    f_ok = compare_filtered(py_dir, jl_dir)
    g_ok = compare_germline(py_dir, jl_dir)

    println("\n══════ Summary ══════")
    println("  filtered.tsv.gz: $(f_ok ? "✓ PASS" : "✗ FAIL")")
    println("  new_V_germline:  $(g_ok ? "✓ PASS" : "✗ FAIL")")

    (f_ok && g_ok) ? println("\n✓ PARITY ACHIEVED") : println("\n✗ PARITY NOT YET ACHIEVED")
    exit(f_ok && g_ok ? 0 : 1)
end

main()
