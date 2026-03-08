#!/usr/bin/env julia
# Run IgDiscover pipeline on ERR1760498 (human Ig heavy chain, paired-end).
#
# Prerequisites: igblastn, makeblastdb, muscle/muscle3, (optional) pear in PATH.
# For faster V-gene discovery, run with threads: julia -t 4 --project=@. scripts/run_err1760498.jl
# Downloads: IMGT human IGH V/D/J, ENA ERR1760498 FASTQ.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using IgDiscover
using FASTX
using CodecZlib

const IMGT_BASE = "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG"
const ENA_BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR176/008/ERR1760498"

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const DB_DIR = joinpath(DATA_DIR, "imgt_human_igh")
const RUN_DIR = joinpath(DATA_DIR, "ERR1760498")
const ANALYSIS_DIR = joinpath(DATA_DIR, "ERR1760498_analysis")

function download_file(url::String, dest::String)
    if isfile(dest)
        @info "Already exists: $dest"
        return
    end
    mkpath(dirname(dest))
    @info "Downloading $url -> $dest"
    cmd = if success(pipeline(`which wget`, devnull))
        `wget -q -O - $url`
    elseif success(pipeline(`which curl`, devnull))
        `curl -sL $url`
    else
        error("Need wget or curl for downloads")
    end
    open(dest, "w") do io
        run(pipeline(cmd, stdout=io))
    end
end

function setup_database()
    mkpath(DB_DIR)
    for (src, dst) in [("IGHV.fasta", "V.fasta"), ("IGHD.fasta", "D.fasta"), ("IGHJ.fasta", "J.fasta")]
        dest = joinpath(DB_DIR, dst)
        if !isfile(dest)
            download_file("$IMGT_BASE/$src", dest)
        end
    end
    @info "Database: $DB_DIR"
end

function fastq_to_fasta(fastq_path::String, fasta_path::String)
    @info "Converting FASTQ -> FASTA: $fastq_path -> $fasta_path"
    mkpath(dirname(fasta_path))
    open(fasta_path, "w") do out
        stream = endswith(fastq_path, ".gz") ?
            FASTQ.Reader(GzipDecompressorStream(open(fastq_path))) :
            open(FASTQ.Reader, fastq_path)
        for rec in stream
            id = first(split(FASTQ.identifier(rec), r"\s+"))
            seq = uppercase(String(FASTQ.sequence(rec)))
            println(out, ">", id)
            println(out, seq)
        end
        close(stream)
    end
    @info "Wrote $fasta_path"
end

function get_reads()
    mkpath(RUN_DIR)
    r1 = joinpath(RUN_DIR, "ERR1760498_1.fastq.gz")
    r2 = joinpath(RUN_DIR, "ERR1760498_2.fastq.gz")
    download_file("$ENA_BASE/ERR1760498_1.fastq.gz", r1)
    download_file("$ENA_BASE/ERR1760498_2.fastq.gz", r2)

    reads_fasta = joinpath(RUN_DIR, "reads.fasta.gz")
    if isfile(reads_fasta)
        @info "Using existing $reads_fasta"
        return reads_fasta
    end

    merged = joinpath(RUN_DIR, "merged.fastq")
    if success(pipeline(`which pear`, devnull))
        nthreads = max(1, Sys.CPU_THREADS)
        @info "Running PEAR to merge paired-end reads (threads=$nthreads)..."
        pear_out = joinpath(RUN_DIR, "pear_out")
        run(pipeline(
            `pear -f $r1 -r $r2 -o $pear_out -j $nthreads`,
            stdout=stdout, stderr=stderr
        ))
        merged_fq = pear_out * ".assembled.fastq"
        isfile(merged_fq) || error("PEAR did not produce $merged_fq")
        fastq_to_fasta(merged_fq, reads_fasta[1:end-3])
        run(`gzip -f $(reads_fasta[1:end-3])`)
    else
        @info "PEAR not found; using R1 only"
        fastq_to_fasta(r1, reads_fasta[1:end-3])
        run(`gzip -f $(reads_fasta[1:end-3])`)
    end
    reads_fasta
end

function main()
    @info "IgDiscover pipeline on ERR1760498"
    setup_database()
    reads = get_reads()
    config_path = joinpath(ANALYSIS_DIR, "igdiscover.toml")
    if isfile(config_path) && isfile(joinpath(ANALYSIS_DIR, "reads.fasta.gz"))
        @info "Analysis dir already initialized, running pipeline..."
    else
        @info "Initializing analysis at $ANALYSIS_DIR"
        IgDiscover.init_analysis(ANALYSIS_DIR, DB_DIR, reads)
    end
    @info "Running pipeline..."
    IgDiscover.run_pipeline(ANALYSIS_DIR)
    @info "Done. Results in $ANALYSIS_DIR"
end

main()
