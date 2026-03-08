# Merge — paired-end read merging via external tools (PEAR or FLASH)
#
# Wraps PEAR or FLASH to merge overlapping paired-end reads into single sequences.
# Supports optional caching of merged files to avoid redundant computation.

# ─── Merge programs as types for dispatch ───

abstract type MergeProgram end

struct PEARMerge <: MergeProgram end
struct FLASHMerge <: MergeProgram
    maximum_overlap::Int
end

FLASHMerge() = FLASHMerge(300)

"""
    ReadMerger(; program=:pear, threads=0, flash_maximum_overlap=300)

Callable struct that merges paired-end reads using an external program.

```julia
merger = ReadMerger(program=:pear, threads=4)
stats = merger("R1.fastq.gz", "R2.fastq.gz", "merged.fastq.gz")
```
"""
struct ReadMerger{P <: MergeProgram}
    program::P
    threads::Int
end

function ReadMerger(; program::Symbol=:pear, threads::Int=0,
                      flash_maximum_overlap::Int=300)
    nthreads = threads > 0 ? threads : Sys.CPU_THREADS
    prog = program === :pear ? PEARMerge() :
           program === :flash ? FLASHMerge(flash_maximum_overlap) :
           error("Unknown merge program: $program. Use :pear or :flash")
    ReadMerger(prog, nthreads)
end

struct MergeStats
    merged::Int
    total::Int
end

# ─── PEAR merging ───

function (m::ReadMerger{PEARMerge})(reads1::AbstractString, reads2::AbstractString,
                                     output::AbstractString)
    Sys.which("pear") !== nothing || error("pear not found in PATH")
    mktempdir() do tmpdir
        prefix = joinpath(tmpdir, "merged")
        cmd = `pear -j $(m.threads) -f $reads1 -r $reads2 -o $prefix`
        log_output = read(cmd, String)
        assembled = prefix * ".assembled.fastq"
        isfile(assembled) || error("PEAR did not produce output")
        compress_fastq(assembled, output)
        parse_pear_stats(log_output)
    end
end

function parse_pear_stats(log::AbstractString)
    merged = 0
    total = 0
    for line in split(log, '\n')
        m_assembled = match(r"Assembled reads \.+: ([\d,]+) / ([\d,]+)", line)
        if m_assembled !== nothing
            merged = parse(Int, replace(m_assembled[1], "," => ""))
            total = parse(Int, replace(m_assembled[2], "," => ""))
        end
    end
    MergeStats(merged, total)
end

# ─── FLASH merging ───

function (m::ReadMerger{FLASHMerge})(reads1::AbstractString, reads2::AbstractString,
                                      output::AbstractString)
    Sys.which("flash") !== nothing || error("flash not found in PATH")
    mktempdir() do tmpdir
        cmd = `flash -M $(m.program.maximum_overlap) -t $(m.threads)
               -d $tmpdir -o merged $reads1 $reads2`
        log_output = read(cmd, String)
        assembled = joinpath(tmpdir, "merged.extendedFrags.fastq")
        isfile(assembled) || error("FLASH did not produce output")
        compress_fastq(assembled, output)
        parse_flash_stats(log_output)
    end
end

function parse_flash_stats(log::AbstractString)
    merged = 0
    total = 0
    for line in split(log, '\n')
        m_combined = match(r"\[FLASH\]\s*Combined pairs:\s*(\d+)", line)
        m_total = match(r"\[FLASH\]\s*Total pairs:\s*(\d+)", line)
        m_combined !== nothing && (merged = parse(Int, m_combined[1]))
        m_total !== nothing && (total = parse(Int, m_total[1]))
    end
    MergeStats(merged, total)
end

# ─── Shared helpers ───

function compress_fastq(input::AbstractString, output::AbstractString)
    gzip = Sys.which("pigz") !== nothing ? "pigz" : "gzip"
    run(`$gzip -c $input` |> output)
end
