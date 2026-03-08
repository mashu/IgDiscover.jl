# Multiple sequence alignment and consensus computation
#
# Aligners are extensible via the Aligner abstract type: define a struct <: Aligner,
# implement run_align(::YourAligner, fasta_str, threads), and register via
# aligner_for_program().

abstract type Aligner end

"""Run the aligner on FASTA string; returns raw alignment output as String."""
function run_align end

struct MafftAligner <: Aligner end
function run_align(::MafftAligner, fasta_str::String, threads::Int)
    read(pipeline(IOBuffer(fasta_str), `mafft --quiet --thread $threads -`), String)
end

struct ClustaloAligner <: Aligner end
function run_align(::ClustaloAligner, fasta_str::String, threads::Int)
    read(pipeline(IOBuffer(fasta_str), `clustalo --threads=$threads --infile=-`), String)
end

# ─── MUSCLE: v3 vs v5 via MuscleVersion dispatch ───
#
# MUSCLE v3 (bioconda default): stdin/stdout, -maxiters, -diags
#   muscle -quiet -maxiters 1 -diags -in /dev/stdin -out /dev/stdout
#
# MUSCLE v5 (separate install):  file-based, -super5 or -align + -refineiters
#   muscle -super5 in.fa -output out.afa -threads N
#
# The two versions have completely incompatible CLIs.

abstract type MuscleVersion end
struct MuscleV3 <: MuscleVersion end
struct MuscleV5 <: MuscleVersion end

struct MuscleAligner{V <: MuscleVersion} <: Aligner
    variant::String   # "muscle", "muscle-fast", "muscle-medium"
end

# Version detection: run once, cache result. For v3, prefer "muscle", fall back to "muscle3".

const DETECTED_MUSCLE_VERSION = Ref{Union{Nothing, MuscleVersion}}(nothing)
const MUSCLE_CMD = Ref{String}("muscle")  # executable name used for alignment (muscle or muscle3)

function detect_muscle_version()
    DETECTED_MUSCLE_VERSION[] !== nothing && return DETECTED_MUSCLE_VERSION[]
    # v3: "MUSCLE v3.8.31 by Robert C. Edgar"
    # v5: "muscle 5.1.xxx" or "muscle 5.2"
    # Use Sys.which to avoid running missing commands; ignorestatus so non-zero exit still returns output.
    if Sys.which("muscle") !== nothing
        output = read(pipeline(ignorestatus(`muscle -version`), stderr=devnull), String)
        if !isempty(output)
            if occursin(r"muscle\s+5\.", output)
                DETECTED_MUSCLE_VERSION[] = MuscleV5()
                MUSCLE_CMD[] = "muscle"
                @info "Detected MUSCLE v5"
            else
                DETECTED_MUSCLE_VERSION[] = MuscleV3()
                MUSCLE_CMD[] = "muscle"
                @info "Detected MUSCLE v3 (muscle)"
            end
            return DETECTED_MUSCLE_VERSION[]
        end
    end
    if Sys.which("muscle3") !== nothing
        output3 = read(pipeline(ignorestatus(`muscle3 -version`), stderr=devnull), String)
        if !isempty(output3)
            DETECTED_MUSCLE_VERSION[] = MuscleV3()
            MUSCLE_CMD[] = "muscle3"
            @info "Detected MUSCLE v3 (muscle3)"
            return DETECTED_MUSCLE_VERSION[]
        end
    end
    error("MUSCLE not found. Install 'muscle' or 'muscle3' (e.g. conda install -c bioconda muscle).")
end

function make_muscle_aligner(variant::String)
    v = detect_muscle_version()
    MuscleAligner{typeof(v)}(variant)
end

# ─── MUSCLE v3: pipe-based, -maxiters, -diags (matches Python exactly) ───

function run_align(a::MuscleAligner{MuscleV3}, fasta_str::String, ::Int)
    exe = MUSCLE_CMD[]
    cmd = if a.variant == "muscle-fast"
        # Python: ["muscle", "-quiet", "-maxiters", "1", "-diags", "-in", "-", "-out", "-"]
        `$exe -quiet -maxiters 1 -diags -in /dev/stdin -out /dev/stdout`
    elseif a.variant == "muscle-medium"
        # Python: ["muscle", "-quiet", "-maxiters", "2", "-diags", "-in", "-", "-out", "-"]
        `$exe -quiet -maxiters 2 -diags -in /dev/stdin -out /dev/stdout`
    else
        # Python: ["muscle", "-quiet", "-in", "-", "-out", "-"]
        `$exe -quiet -in /dev/stdin -out /dev/stdout`
    end
    read(pipeline(IOBuffer(fasta_str), cmd), String)
end

# ─── MUSCLE v5: file-based, -super5 for fast, -align + -refineiters for others ───
#
# v5 has no -maxiters or -diags.
# -super5 is the fast mode (HMM ensemble, closest to v3's -maxiters 1 -diags).
# -align with -refineiters 0 gives progressive alignment only (no refinement).

function run_align(a::MuscleAligner{MuscleV5}, fasta_str::String, threads::Int)
    tmpin = tempname() * ".fa"
    tmpout = tempname() * ".afa"
    write(tmpin, fasta_str)

    cmd = if a.variant == "muscle-fast"
        # -super5: fastest v5 mode, HMM-based, no iterative refinement
        `muscle -super5 $tmpin -output $tmpout -threads $threads`
    elseif a.variant == "muscle-medium"
        # Progressive alignment + 1 refinement pass
        `muscle -align $tmpin -output $tmpout -refineiters 1 -threads $threads`
    else
        # Full: progressive + many refinement passes
        `muscle -align $tmpin -output $tmpout -refineiters 100 -threads $threads`
    end

    run(pipeline(cmd, stderr=devnull))
    result = read(tmpout, String)
    rm(tmpin; force=true)
    rm(tmpout; force=true)
    result
end

# ─── Aligner registry (lazy — version detection at first use) ───

const ALIGNERS = Dict{String, Aligner}()

function aligner_for_program(program::String)
    get!(ALIGNERS, program) do
        if program in ("muscle", "muscle-fast", "muscle-medium")
            make_muscle_aligner(program)
        elseif program == "mafft"
            MafftAligner()
        elseif program == "clustalo"
            ClustaloAligner()
        else
            error("Alignment program '$program' not supported. " *
                  "Known: mafft, clustalo, muscle, muscle-fast, muscle-medium")
        end
    end
end

"""
    parse_fasta_string(output) -> Dict{String,String}

Parse a FASTA-format string into name → sequence dictionary.
"""
parse_fasta_string(output::String) = read_fasta_dict_from_io(IOBuffer(output))

"""
    multialign(sequences; program="muscle-fast", threads=Sys.CPU_THREADS) -> Dict{String,String}

Run a multiple sequence alignment. Dispatch is via the Aligner abstraction.
"""
function multialign(sequences::Dict{String,String};
                    program::String = "muscle-fast",
                    threads::Int = Sys.CPU_THREADS)
    isempty(sequences) && return Dict{String,String}()
    fasta_input = IOBuffer()
    for (name, seq) in sequences
        println(fasta_input, ">", name)
        println(fasta_input, seq)
    end
    fasta_str = String(take!(fasta_input))
    output = run_align(aligner_for_program(program), fasta_str, threads)
    parse_fasta_string(output)
end

# ─── Nucleotide counter (replaces Dict{Char,Int} in consensus hot loop) ───

mutable struct NucleotideCounter
    a::Int; c::Int; g::Int; t::Int; n::Int; gap::Int
end

NucleotideCounter() = NucleotideCounter(0, 0, 0, 0, 0, 0)

function reset!(nc::NucleotideCounter)
    nc.a = nc.c = nc.g = nc.t = nc.n = nc.gap = 0
    nc
end

function count!(nc::NucleotideCounter, ch::Char)
    ch == 'A' ? (nc.a += 1) :
    ch == 'C' ? (nc.c += 1) :
    ch == 'G' ? (nc.g += 1) :
    ch == 'T' ? (nc.t += 1) :
    ch == '-' ? (nc.gap += 1) :
    (nc.n += 1)
    nc
end

function best_base(nc::NucleotideCounter)
    best_char, best_freq = '-', nc.gap
    nc.a > best_freq && (best_char = 'A'; best_freq = nc.a)
    nc.c > best_freq && (best_char = 'C'; best_freq = nc.c)
    nc.g > best_freq && (best_char = 'G'; best_freq = nc.g)
    nc.t > best_freq && (best_char = 'T'; best_freq = nc.t)
    nc.n > best_freq && (best_char = 'N'; best_freq = nc.n)
    (best_char, best_freq)
end

function adjust_gaps!(nc::NucleotideCounter, excess::Int)
    nc.gap = excess > 0 ? excess : 0
    nc
end

"""
    consensus_sequence(aligned; threshold=0.7, ambiguous='N') -> String

Compute a consensus from aligned sequences, allowing degraded 3' ends.
"""
function consensus_sequence(aligned::AbstractVector{String};
                           threshold::Float64 = 0.7,
                           ambiguous::Char = 'N')
    isempty(aligned) && return ""
    ncols = maximum(length, aligned)
    nseqs = length(aligned)

    result_rev = Char[]
    sizehint!(result_rev, ncols)
    active = max(1, round(Int, nseqs * 0.05))

    counts = NucleotideCounter()

    for col_rev in 0:ncols-1
        col = ncols - col_rev
        reset!(counts)

        for seq in aligned
            ch = col <= length(seq) ? seq[col] : '-'
            count!(counts, ch)
        end

        gap_count = counts.gap
        active = max(nseqs - gap_count, active)

        adjust_gaps!(counts, gap_count - (nseqs - active))

        if col_rev >= 10
            active = nseqs
        end

        ch, freq = best_base(counts)

        if freq / active >= threshold
            ch != '-' && push!(result_rev, ch)
        else
            push!(result_rev, ambiguous)
        end
    end
    String(reverse!(result_rev))
end

function consensus_sequence(aligned::Dict{String,String}; kwargs...)
    consensus_sequence(collect(values(aligned)); kwargs...)
end

"""
    iterative_consensus(sequences; program, threshold, subsample_size, maximum_subsample_size) -> String

Compute consensus by iteratively increasing subsample size until no N bases remain.
"""
function iterative_consensus(sequences::Vector{String};
                             program::String = "muscle-medium",
                             threshold::Float64 = 0.6,
                             subsample_size::Int = 200,
                             maximum_subsample_size::Int = 1600,
                             rng::Random.AbstractRNG = Random.GLOBAL_RNG)
    while true
        sample = reservoir_sample(sequences, subsample_size; rng = rng)
        seqs = Dict(string(i) => s for (i, s) in enumerate(sample))
        aligned = multialign(seqs; program = program)
        cons = strip(consensus_sequence(aligned; threshold = threshold), 'N')

        !occursin('N', cons) && return cons
        length(sequences) <= subsample_size && return cons
        subsample_size *= 2
        subsample_size > maximum_subsample_size && return cons
    end
end

"""
    reservoir_sample(population, k; rng=Random.GLOBAL_RNG) -> Vector

Reservoir sampling: return at most k randomly chosen elements from population.
"""
function reservoir_sample(population::AbstractVector{T}, k::Int;
                         rng::Random.AbstractRNG = Random.GLOBAL_RNG) where T
    n = length(population)
    k >= n && return copy(population)
    sample = population[1:k]
    for i in (k+1):n
        j = rand(rng, 1:i)
        if j <= k
            sample[j] = population[i]
        end
    end
    sample
end

# ─── Pairwise alignment via BioAlignments ───

struct Alignment
    ref_row::String
    query_row::String
    score::Int
    errors::Int
end

"""
    align_affine(ref, query; gap_open=-6, gap_extend=-1, mismatch=-3, match_score=1) -> Alignment

Global alignment with affine gap penalties (BWA-like defaults).
The gap penalty convention is: gap of length k costs `gap_open + (k-1)*gap_extend`.
BioAlignments uses `bio_gap_open + k*gap_extend`, so we adjust accordingly.
"""
function align_affine(ref::AbstractString, query::AbstractString;
                     gap_open::Int = -6, gap_extend::Int = -1,
                     mismatch::Int = -3, match_score::Int = 1)
    ref_dna = BioSequences.LongDNA{4}(ref)
    query_dna = BioSequences.LongDNA{4}(query)

    bio_gap_open = gap_open - gap_extend
    model = BioAlignments.AffineGapScoreModel(
        match=match_score, mismatch=mismatch,
        gap_open=bio_gap_open, gap_extend=gap_extend,
    )
    result = BioAlignments.pairalign(BioAlignments.GlobalAlignment(), query_dna, ref_dna, model)
    aln = BioAlignments.alignment(result)

    ref_chars = Char[]
    query_chars = Char[]
    errors = 0
    for (seq_nuc, ref_nuc) in aln
        sc = Char(seq_nuc)
        rc = Char(ref_nuc)
        push!(ref_chars, rc)
        push!(query_chars, sc)
        sc != rc && (errors += 1)
    end

    Alignment(String(ref_chars), String(query_chars),
              BioAlignments.score(result), errors)
end

# ─── describe_nt_change: column kinds dispatched ───

abstract type ColumnKind end
struct MatchCol <: ColumnKind end
struct InsertCol <: ColumnKind end
struct DeleteCol <: ColumnKind end
struct SubstCol <: ColumnKind end

column_kind(c1::Char, c2::Char) =
    c1 == c2 ? MatchCol() : (c1 == '-' ? InsertCol() : (c2 == '-' ? DeleteCol() : SubstCol()))

"""Process one column (or run for indels); return (next_idx, next_i)."""
function apply_column! end

apply_column!(changes::Vector{String}, idx::Int, i::Int, aln::Alignment, ::MatchCol) =
    (idx + 1, i + 1)

function apply_column!(changes::Vector{String}, idx::Int, i::Int, aln::Alignment, ::InsertCol)
    buf = IOBuffer()
    i_end = i
    while i_end <= length(aln.ref_row) && aln.ref_row[i_end] == '-'
        write(buf, aln.query_row[i_end])
        i_end += 1
    end
    push!(changes, "$(idx-1)_$(idx)ins$(String(take!(buf)))")
    (idx, i_end)
end

function apply_column!(changes::Vector{String}, idx::Int, i::Int, aln::Alignment, ::DeleteCol)
    buf = IOBuffer()
    start_idx = idx
    i_end = i
    while i_end <= length(aln.ref_row) && aln.query_row[i_end] == '-'
        write(buf, aln.ref_row[i_end])
        i_end += 1
        idx += 1
    end
    push!(changes, "$(start_idx)_$(idx-1)del$(String(take!(buf)))")
    (idx, i_end)
end

apply_column!(changes::Vector{String}, idx::Int, i::Int, aln::Alignment, ::SubstCol) = begin
    c1, c2 = aln.ref_row[i], aln.query_row[i]
    push!(changes, "$(idx)$(c1)>$(c2)")
    (idx + 1, i + 1)
end

"""
    describe_nt_change(ref, query) -> String

Describe nucleotide changes between ref and query using HGVS-like notation.
Column type (match/insert/delete/subst) is dispatched via ColumnKind.
"""
function describe_nt_change(ref::AbstractString, query::AbstractString)
    aln = align_affine(ref, query)
    changes = String[]
    idx = 1
    i = 1
    while i <= length(aln.ref_row)
        kind = column_kind(aln.ref_row[i], aln.query_row[i])
        idx, i = apply_column!(changes, idx, i, aln, kind)
    end
    join(changes, "; ")
end
