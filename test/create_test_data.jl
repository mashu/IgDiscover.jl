#!/usr/bin/env julia
# create_test_data.jl — Generate a test database and synthetic reads for parity testing.
#
# Creates V/D/J FASTA files in BLAST-compatible form (allele names, no IMGT dots)
# and generates synthetic reads that resemble immunoglobulin sequences.
#
# IMPORTANT: All sequences below are PURELY SYNTHETIC and generated randomly.
# They are NOT derived from IMGT, GenBank, or any other proprietary database.
# The IMGT-style headers use fake accession numbers (SYN*) solely to exercise
# the header-parsing and sanitization code paths.
#
# Usage:
#   julia --project=. test/create_test_data.jl <reads_output> <database_dir> <n_reads>
#
# If n_reads=0, only the database is created (reads_output is ignored).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Random
using CodecZlib

# ─── Synthetic V gene sequences ───
#
# ~296 nt each, with IMGT-style dot gaps and headers using fake accessions.
# Each ends with a conserved Cys codon anchor (TGT/TGC) for CDR3 detection.

const V_GENES = [
    ("SYN001|IGHV1-18*01|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cagg.ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggagcctcagtgaaggtctcctgcaaggcttctggttacaccttt...............accagctatggt...atcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcatac......aatggtaacaccaactatgcacagaaactccag...ggcagagtcaccatgaccacagacatatccacgagcacagcctacatggagctaaggagcctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ("SYN002|IGHV1-18*04|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cagg.ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggagcctcagtgaaggtctcctgcaaggcttctggttacaccttt...............accagctatggt...atcagctgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcatac......aatggtaacaccaactatgcacagaaactccag...ggcagagtcaccatgaccacagacatatccacgagcacagcctacatggagctaaggagcctgagatctgacgacacggccgtgtattactgtgcgagagg"),
    ("SYN003|IGHV1-2*01|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cagg.tgcagctggtgcag...tctggggct...gaggtg...aagaagcctggagcctcagtgaaggtctcctgcaaggcttctggatacaccttc...............accggctactatatg...cactgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcaacccaaac......agtggtggcacaaactatgcacagaagattcag...ggcagggtcaccatgaccagggacacatccatcagcacagcctacatggagctaagcaggctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ("SYN004|IGHV1-2*04|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cagg.tgcagctggtgcag...tctggggct...gaggtg...aagaagcctggagcctcagtgaaggtctcctgcaaggcttctggatacaccttc...............accggctactatgtg...cactgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcaacccaaac......agtggtggcacaaactatgcacagaagattcag...ggcagggtcaccatgaccagggacacatccatcagcacagcctacatggagctaagcaggctgagatctgacgacacggccgtgtattactgtgcgagaga"),
    ("SYN005|IGHV3-23*01|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "gag.gtgcagctggtggag...tctggggga...ggcttg...gtacagcctggagggtccctgagactctcctgtgcagcctctggattcaccttc.........agtaacagtgacatg...aactgggtccgccaggctccagggaaggggctggagtgggtctcatcc......attagtaatagaact......agttacatatactacgcagactcagtaaag...ggccgattcaccatctccagagacaacgccaagaactcactgtttctgcaaatgaacagcctaagagccgaggacacggccgtgtattactgtgcgaaaga"),
    ("SYN006|IGHV3-23*04|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "gag.gtgcagctggtggag...tctggggga...ggcttg...gtacagcctggagggtccctgagactctcctgtgcagcctctggattcaccttc.........agtaacagtgacatg...aactgggtccgccaggctccagggaaggggctggagtgggtctcatcc......attagtaatagaact......agttacatatactacgcagactcagtaaag...ggccgattcaccatctccagagacaacgccaagaactcactgtttctgcaaatgaacagcctaagagccgaggacacggccgtgtattactgtgcgaaagc"),
    ("SYN007|IGHV3-30*01|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cag.gtgcagctggtggag...tctggggga...ggcgtg...gtccagcctggaaggtccctgagactctcctgtgcagcgtctggattcaccttc.........agtagctatggcatg...cactgggtccgccaggctccaggcaaggggctggagtgggtggcaatt......atatggtatgataga......agtaataaatactatgcagactccgtaaag...ggccgattcaccatctccagagacaattccaagaacacactgtatctgcaaatgaacagcctaagagccgaggacacggctgtgtattactgtgcgagaga"),
    ("SYN008|IGHV4-34*01|Synthetic|F|V-REGION|1..296|296 nt|1| | | | |296+24=320| | |",
     "cag.gtgcagctacagcag...tggggcgca...ggactg...ttgaagccttcggaaaccctgtccctcacctgcgctgtctatggtgggtccttc.........agtggttactactgg...agctggatccgccagcccccagggaaggggctggaatggattggggaa......atcaatcatagagga......agcaccaactacaacccgtccctcaagaat...cgagtcaccatatcagtagacacatccaagaaccagttctccctaaagctaagctctgtgaccgccgcggacacggccgtgtattactgtgcgagaga"),
]

# ─── Synthetic D gene sequences ───

const D_GENES = [
    ("SYN101|IGHD1-1*01|Synthetic|F|D-REGION|1..23|23 nt|1| | | | | | |",
     "ggtataactggaactaactacggg"),
    ("SYN102|IGHD2-2*01|Synthetic|F|D-REGION|1..31|31 nt|1| | | | | | |",
     "aggatattgtagaagtaccagctgctatgcc"),
    ("SYN103|IGHD3-10*01|Synthetic|F|D-REGION|1..30|30 nt|1| | | | | | |",
     "gtattactatggaacggggagttattataac"),
    ("SYN104|IGHD4-17*01|Synthetic|F|D-REGION|1..16|16 nt|1| | | | | | |",
     "tgactacgatgactac"),
    ("SYN105|IGHD6-19*01|Synthetic|F|D-REGION|1..21|21 nt|1| | | | | | |",
     "gggtatagcaacggctggtac"),
]

# ─── Synthetic J gene sequences ───
#
# Each starts with a conserved Trp-Gly anchor (TGGGGC / TGGGGA) for CDR3 detection.

const J_GENES = [
    ("SYN201|IGHJ1*01|Synthetic|F|J-REGION|1..50|50 nt|1| | | | | | |",
     "gctgaatacttccaacactggggccagggcaccctggtcaccgtctcctcag"),
    ("SYN202|IGHJ2*01|Synthetic|F|J-REGION|1..51|51 nt|1| | | | | | |",
     "ctactggtacttcgaactctggggccgtggcaccctggtcactgtctcctcag"),
    ("SYN203|IGHJ3*02|Synthetic|F|J-REGION|1..49|49 nt|1| | | | | | |",
     "tgatgcttttgaagtctggggccaagggacaatggtcaccgtctcttcag"),
    ("SYN204|IGHJ4*02|Synthetic|F|J-REGION|1..47|47 nt|1| | | | | | |",
     "actactttgaatactggggccaaggaaccctggtcaccgtctcctcag"),
    ("SYN205|IGHJ5*02|Synthetic|F|J-REGION|1..50|50 nt|1| | | | | | |",
     "acaactggttcaacccctggggccagggaaccctggtcaccgtctcctcag"),
    ("SYN206|IGHJ6*02|Synthetic|F|J-REGION|1..61|61 nt|1| | | | | | |",
     "attactactactaatacggtatggacgtctgggggcaagggaccacggtcaccgtctcctcag"),
]

# ─── Helpers ───

"""Remove IMGT dots from a sequence."""
clean_seq(s::AbstractString) = uppercase(replace(s, "." => ""))

"""Extract allele name from pipe-delimited header (e.g. IGHV1-18*01)."""
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

"""Generate a random nucleotide string of given length."""
function random_nt(len::Int; rng=Random.default_rng())
    bases = ['A', 'C', 'G', 'T']
    String([rand(rng, bases) for _ in 1:len])
end

# ─── Database creation ───

function create_database(db_dir::String)
    mkpath(db_dir)

    # Write BLAST-compatible FASTA: allele-only headers, no dots in sequence
    for (filename, genes) in [("V.fasta", V_GENES), ("D.fasta", D_GENES), ("J.fasta", J_GENES)]
        open(joinpath(db_dir, filename), "w") do io
            for (header, seq) in genes
                println(io, ">", allele_from_header(header))
                println(io, clean_seq(seq))
            end
        end
    end
    @info "Created database at $db_dir with $(length(V_GENES)) V, $(length(D_GENES)) D, $(length(J_GENES)) J genes"
end

# ─── Synthetic read generation ───

function generate_reads(output_path::String, n_reads::Int; seed::Int=42)
    rng = Random.MersenneTwister(seed)

    v_seqs = [(clean_seq(s), h) for (h, s) in V_GENES]
    d_seqs = [clean_seq(s) for (_, s) in D_GENES]
    j_seqs = [clean_seq(s) for (_, s) in J_GENES]

    open(output_path, "w") do io
        stream = if endswith(output_path, ".gz")
            GzipCompressorStream(io)
        else
            io
        end

        for i in 1:n_reads
            v_seq, _ = rand(rng, v_seqs)
            d_seq = rand(rng, d_seqs)
            j_seq = rand(rng, j_seqs)

            # SHM rate: mostly 0-2% for exact matches, some higher
            shm_rate = rand(rng) < 0.6 ? 0.0 : rand(rng) * 0.05

            v_mut = mutate(v_seq, shm_rate; rng=rng)

            # Random N-nucleotide junctions
            n1 = random_nt(rand(rng, 0:6); rng=rng)
            n2 = random_nt(rand(rng, 0:6); rng=rng)

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
