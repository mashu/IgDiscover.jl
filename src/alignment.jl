# Multiple sequence alignment and consensus computation
#
# Aligners are extensible via the Aligner abstract type: define a struct <: Aligner,
# implement run_align(::YourAligner, fasta_str, threads), and register in ALIGNERS.

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

struct MuscleAligner <: Aligner
    variant::String  # "muscle", "muscle-fast", "muscle-medium"
end
function run_align(a::MuscleAligner, fasta_str::String, ::Int)
    tmpin = tempname() * ".fa"
    tmpout = tempname() * ".afa"
    write(tmpin, fasta_str)
    p = run(pipeline(`muscle -align $tmpin -output $tmpout`, stderr=devnull); wait=false)
    wait(p)
    result = if p.exitcode == 0 && isfile(tmpout) && filesize(tmpout) > 0
        read(tmpout, String)
    else
        rm(tmpout; force=true)
        run(pipeline(`muscle -quiet -in $tmpin -out $tmpout`, stderr=devnull))
        read(tmpout, String)
    end
    rm(tmpin; force=true)
    rm(tmpout; force=true)
    result
end

const ALIGNERS = Dict{String, Aligner}(
    "mafft" => MafftAligner(),
    "clustalo" => ClustaloAligner(),
    "muscle" => MuscleAligner("muscle"),
    "muscle-fast" => MuscleAligner("muscle-fast"),
    "muscle-medium" => MuscleAligner("muscle-medium"),
)

aligner_for_program(program::String) = get(ALIGNERS, program) do
    error("Alignment program '$program' not supported. Registered: ", join(sort!(collect(keys(ALIGNERS))), ", "))
end

"""
    parse_fasta_string(output) -> Dict{String,String}

Parse a FASTA-format string into name → sequence dictionary. Delegates to read_fasta_dict_from_io.
"""
parse_fasta_string(output::String) = read_fasta_dict_from_io(IOBuffer(output))

"""
    multialign(sequences::Dict{String,String}; program="muscle-fast", threads=Sys.CPU_THREADS) -> Dict{String,String}

Run a multiple sequence alignment program and return aligned sequences.
Dispatch is via the `Aligner` abstraction; register new programs in `ALIGNERS` and implement `run_align`.
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
    if excess > 0
        nc.gap = excess
    else
        nc.gap = 0
    end
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
                             maximum_subsample_size::Int = 1600)
    while true
        sample = reservoir_sample(sequences, subsample_size)
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
    reservoir_sample(population, k) -> Vector

Reservoir sampling: return at most k randomly chosen elements from population.
"""
function reservoir_sample(population::AbstractVector{T}, k::Int) where T
    n = length(population)
    k >= n && return copy(population)
    sample = population[1:k]
    for i in (k+1):n
        j = rand(1:i)
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
Delegates to BioAlignments.pairalign; returns aligned strings with gaps for downstream use.

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
