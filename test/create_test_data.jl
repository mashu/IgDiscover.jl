#!/usr/bin/env julia
# create_test_data.jl — Generate a test database and synthetic reads for parity testing.
#
# Creates V/D/J FASTA files in BLAST-compatible form (allele names, no IMGT dots)
# and generates synthetic reads that resemble immunoglobulin sequences.
#
# Usage:
#   julia --project=. test/create_test_data.jl <reads_output> <database_dir> <n_reads>
#
# If n_reads=0, only the database is created (reads_output is ignored).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using CodecZlib

# ─── IMGT-style V gene sequences (realistic human IGHV, with dots) ───

const V_GENES = [
    ("M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |",
     "cagg.ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggttacaccttt...............accagctatggt...atcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcttac......aatggtaacacaaactatgcacagaagctccag...ggcagagtcaccatgaccacagacacatccacgagcacagcctacatggagctgaggagcctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ("M99641|IGHV1-18*04|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |",
     "cagg.ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggttacaccttt...............accagctatggt...atcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcttac......aatggtaacacaaactatgcacagaagctccag...ggcagagtcaccatgaccacagacacatccacgagcacagcctacatggagctgaggagcctgagatctgacgacacggccgtgtattactgtgcgagagg"),
    ("X62106|IGHV1-2*01|Homo sapiens|F|V-REGION|293..588|296 nt|1| | | | |296+24=320| | |",
     "cagg.tgcagctggtgcag...tctggggct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggatacaccttc...............accggctactatatg...cactgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcaaccctaac......agtggtggcacaaactatgcacagaagtttcag...ggcagggtcaccatgaccagggacacgtccatcagcacagcctacatggagctgagcaggctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ("AB019438|IGHV1-2*04|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cagg.tgcagctggtgcag...tctggggct...gaggtg...aagaagcctggggcctcagtgaaggtctcctgcaaggcttctggatacaccttc...............accggctactatgtg...cactgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcaaccctaac......agtggtggcacaaactatgcacagaagtttcag...ggcagggtcaccatgaccagggacacgtccatcagcacagcctacatggagctgagcaggctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ("M99649|IGHV3-23*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "gag.gtgcagctggtggag...tctggggga...ggcttg...gtacagcctggggggtccctgagactctcctgtgcagcctctggattcaccttc.........agtaacagtgacatg...aactgggtccgccaggctccagggaaggggctggagtgggtctcatcc......attagtagtagtact......agttacatatactacgcagactcagtgaag...ggccgattcaccatctccagagacaacgccaagaactcactgtatctgcaaatgaacagcctgagagccgaggacacggccgtgtattactgtgcgaaaga"),
    ("M99649|IGHV3-23*04|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "gag.gtgcagctggtggag...tctggggga...ggcttg...gtacagcctggggggtccctgagactctcctgtgcagcctctggattcaccttc.........agtaacagtgacatg...aactgggtccgccaggctccagggaaggggctggagtgggtctcatcc......attagtagtagtact......agttacatatactacgcagactcagtgaag...ggccgattcaccatctccagagacaacgccaagaactcactgtatctgcaaatgaacagcctgagagccgaggacacggccgtgtattactgtgcgaaagc"),
    ("AB019439|IGHV3-30*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cag.gtgcagctggtggag...tctggggga...ggcgtg...gtccagcctgggaggtccctgagactctcctgtgcagcgtctggattcaccttc.........agtagctatggcatg...cactgggtccgccaggctccaggcaaggggctggagtgggtggcagtt......atatggtatgatgga......agtaataaatactatgcagactccgtgaag...ggccgattcaccatctccagagacaattccaagaacacgctgtatctgcaaatgaacagcctgagagccgaggacacggctgtgtattactgtgcgagaga"),
    ("X92218|IGHV4-34*01|Homo sapiens|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cag.gtgcagctacagcag...tggggcgca...ggactg...ttgaagccttcggagaccctgtccctcacctgcgctgtctatggtgggtccttc.........agtggttactactgg...agctggatccgccagcccccagggaaggggctggagtggattggggaa......atcaatcatagtgga......agcaccaactacaacccgtccctcaagagt...cgagtcaccatatcagtagacacgtccaagaaccagttctccctgaagctgagctctgtgaccgccgcggacacggccgtgtattactgtgcgagaga"),
]

const D_GENES = [
    ("X97051|IGHD1-1*01|Homo sapiens|F|D-REGION|8..30|23 nt|1| | | | | | |",
     "ggtataactggaactgactacggg"),
    ("X13972|IGHD2-2*01|Homo sapiens|F|D-REGION|8..25|18 nt|1| | | | | | |",
     "aggatattgtagtagtaccagctgctatgcc"),
    ("X97051|IGHD3-10*01|Homo sapiens|F|D-REGION|8..22|15 nt|1| | | | | | |",
     "gtattactatggttcggggagttattataac"),
    ("X97051|IGHD4-17*01|Homo sapiens|F|D-REGION|8..22|15 nt|1| | | | | | |",
     "tgactacggtgactac"),
    ("X97051|IGHD6-19*01|Homo sapiens|F|D-REGION|8..24|17 nt|1| | | | | | |",
     "gggtatagcagcggctggtac"),
]

const J_GENES = [
    ("J00256|IGHJ1*01|Homo sapiens|F|J-REGION|9..53|45 nt|1| | | | | | |",
     "gctgaatacttccagcactggggccagggcaccctggtcaccgtctcctcag"),
    ("J00256|IGHJ2*01|Homo sapiens|F|J-REGION|118..170|53 nt|1| | | | | | |",
     "ctactggtacttcgatctctggggccgtggcaccctggtcactgtctcctcag"),
    ("J00256|IGHJ3*02|Homo sapiens|F|J-REGION|229..279|51 nt|1| | | | | | |",
     "tgatgcttttgatgtctggggccaagggacaatggtcaccgtctcttcag"),
    ("J00256|IGHJ4*02|Homo sapiens|F|J-REGION|346..393|48 nt|1| | | | | | |",
     "actactttgactactggggccaaggaaccctggtcaccgtctcctcag"),
    ("J00256|IGHJ5*02|Homo sapiens|F|J-REGION|459..507|49 nt|1| | | | | | |",
     "acaactggttcgacccctggggccagggaaccctggtcaccgtctcctcag"),
    ("J00256|IGHJ6*02|Homo sapiens|F|J-REGION|572..625|54 nt|1| | | | | | |",
     "attactactactactacggtatggacgtctgggggcaagggaccacggtcaccgtctcctcag"),
]

# ─── Helpers ───

"""Remove IMGT dots from a sequence."""
clean_seq(s::AbstractString) = uppercase(replace(s, "." => ""))

"""Extract allele name from IMGT pipe-delimited header (e.g. IGHV1-18*01) for BLAST-safe FASTA."""
allele_from_header(h::AbstractString) = (p = split(h, '|'); length(p) >= 2 ? String(p[2]) : String(h))

"""Introduce random point mutations at a given rate."""
function mutate(seq::String, rate::Float64; rng=Random.default_rng())
    bases = ['A', 'C', 'G', 'T']
    chars = collect(seq)
    for i in eachindex(chars)
        if rand(rng) < rate
            chars[i] = rand(rng, filter(!=(chars[i]), bases))
        end
    end
    String(chars)
end

"""Generate a random CDR3 of given nucleotide length (must be divisible by 3)."""
function random_cdr3(len::Int; rng=Random.default_rng())
    bases = ['A', 'C', 'G', 'T']
    String([rand(rng, bases) for _ in 1:len])
end

# ─── Database creation ───

function create_database(db_dir::String)
    mkpath(db_dir)

    # Write BLAST-compatible FASTA: allele-only headers, no IMGT dots in sequence
    # (makeblastdb rejects dots and pipe-heavy headers)
    for (filename, genes) in [("V.fasta", V_GENES), ("D.fasta", D_GENES), ("J.fasta", J_GENES)]
        open(joinpath(db_dir, filename), "w") do io
            for (header, seq) in genes
                println(io, ">", allele_from_header(header))
                println(io, clean_seq(seq))
            end
        end
    end
    @info "Created BLAST-compatible database at $db_dir with $(length(V_GENES)) V, $(length(D_GENES)) D, $(length(J_GENES)) J genes"
end

# ─── Synthetic read generation ───

function generate_reads(output_path::String, n_reads::Int; seed::Int=42)
    rng = Random.MersenneTwister(seed)

    # Clean sequences for read generation
    v_seqs = [(clean_seq(s), h) for (h, s) in V_GENES]
    d_seqs = [clean_seq(s) for (_, s) in D_GENES]
    j_seqs = [clean_seq(s) for (_, s) in J_GENES]

    open(output_path, "w") do io
        # Use GzipCompressorStream if output is .gz
        stream = if endswith(output_path, ".gz")
            GzipCompressorStream(io)
        else
            io
        end

        for i in 1:n_reads
            # Pick random V, D, J
            v_seq, _ = rand(rng, v_seqs)
            d_seq = rand(rng, d_seqs)
            j_seq = rand(rng, j_seqs)

            # SHM rate: mostly 0-2% for exact matches, some higher
            shm_rate = rand(rng) < 0.6 ? 0.0 : rand(rng) * 0.05

            # Mutate V region
            v_mut = mutate(v_seq, shm_rate; rng=rng)

            # Random N-nucleotide junctions
            n1 = random_cdr3(rand(rng, 0:6); rng=rng)
            n2 = random_cdr3(rand(rng, 0:6); rng=rng)

            # D region (sometimes trimmed)
            d_trim_5 = rand(rng, 0:min(5, length(d_seq) ÷ 3))
            d_trim_3 = rand(rng, 0:min(5, length(d_seq) ÷ 3))
            d_used = d_seq[1+d_trim_5:end-d_trim_3]

            # Assemble full read: V + N1 + D + N2 + J
            read_seq = v_mut * n1 * d_used * n2 * j_seq

            println(stream, ">read_$(lpad(i, 6, '0'))")
            println(stream, read_seq)
        end

        if stream !== io
            close(stream)
        end
    end
    @info "Generated $n_reads synthetic reads to $output_path"
end

# ─── Main ───

function main()
    if length(ARGS) < 3
        println("Usage: julia create_test_data.jl <reads_output> <database_dir> <n_reads>")
        println("  If n_reads=0, only database is created.")
        exit(1)
    end

    reads_output = ARGS[1]
    db_dir = ARGS[2]
    n_reads = parse(Int, ARGS[3])

    create_database(db_dir)

    if n_reads > 0
        generate_reads(reads_output, n_reads)
    end
end

main()
