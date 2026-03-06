# Multiple sequence alignment and consensus computation

"""
    multialign(sequences::Dict{String,String}; program="muscle-fast", threads=Sys.CPU_THREADS) -> Dict{String,String}

Run a multiple sequence alignment program and return aligned sequences.
Supported programs: "muscle", "muscle-fast", "muscle-medium", "mafft", "clustalo".
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

    output = if program == "mafft"
        read(pipeline(IOBuffer(fasta_str), `mafft --quiet --thread $threads -`), String)
    elseif program == "clustalo"
        read(pipeline(IOBuffer(fasta_str), `clustalo --threads=$threads --infile=-`), String)
    elseif program in ("muscle", "muscle-fast", "muscle-medium")
        tmpin = tempname() * ".fa"
        tmpout = tempname() * ".afa"
        write(tmpin, fasta_str)
        # Try MUSCLE 5 syntax first
        p = run(pipeline(`muscle -align $tmpin -output $tmpout`, stderr=devnull); wait=false)
        wait(p)
        result = if p.exitcode == 0 && isfile(tmpout) && filesize(tmpout) > 0
            read(tmpout, String)
        else
            # Clean stale output before MUSCLE 3 fallback
            rm(tmpout; force=true)
            run(pipeline(`muscle -quiet -in $tmpin -out $tmpout`, stderr=devnull))
            read(tmpout, String)
        end
        rm(tmpin; force=true)
        rm(tmpout; force=true)
        result
    else
        error("Alignment program '$program' not supported")
    end

    parse_fasta_string(output)
end

function parse_fasta_string(output::String)
    aligned = Dict{String,String}()
    current_name = ""
    current_seq = IOBuffer()
    for line in eachline(IOBuffer(output))
        if startswith(line, '>')
            if !isempty(current_name)
                aligned[current_name] = uppercase(String(take!(current_seq)))
            end
            current_name = strip(line[2:end])
        else
            write(current_seq, strip(line))
        end
    end
    if !isempty(current_name)
        aligned[current_name] = uppercase(String(take!(current_seq)))
    end
    aligned
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

    # Build result in reverse order (3'→5'), then reverse at end — O(n) vs O(n²)
    result_rev = Char[]
    sizehint!(result_rev, ncols)
    active = max(1, round(Int, nseqs * 0.05))

    counts = Dict{Char,Int}()

    # Process from the 3' end (reversed), same as Python version
    for col_rev in 0:ncols-1
        col = ncols - col_rev
        empty!(counts)
        for seq in aligned
            c = col <= length(seq) ? seq[col] : '-'
            counts[c] = get(counts, c, 0) + 1
        end
        gap_count = get(counts, '-', 0)
        active = max(nseqs - gap_count, active)

        adjusted_gaps = gap_count - (nseqs - active)
        if adjusted_gaps > 0
            counts['-'] = adjusted_gaps
        else
            delete!(counts, '-')
        end

        if col_rev >= 10
            active = nseqs
        end

        best_char, best_freq = ' ', 0
        for (c, f) in counts
            if f > best_freq
                best_char, best_freq = c, f
            end
        end

        if best_freq / active >= threshold
            if best_char != '-'
                push!(result_rev, best_char)
            end
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
    iterative_consensus(sequences::Vector{String};
                       program="muscle-medium", threshold=0.6,
                       subsample_size=200, maximum_subsample_size=1600) -> String

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

# ─── Affine gap alignment for describe_nt_change ───

struct Alignment
    ref_row::String
    query_row::String
    score::Int
    errors::Int
end

"""
    align_affine(ref, query; gap_open=-6, gap_extend=-1, mismatch=-3, match_score=1) -> Alignment

Global alignment with affine gap penalties (BWA defaults).
"""
function align_affine(ref::AbstractString, query::AbstractString;
                     gap_open::Int = -6, gap_extend::Int = -1,
                     mismatch::Int = -3, match_score::Int = 1)
    m, n = length(ref), length(query)
    INF = typemin(Int) ÷ 2

    # Matrices: M=match/mismatch, H=horizontal gap (in ref), V=vertical gap (in query)
    M = zeros(Int, m + 1, n + 1)
    H = zeros(Int, m + 1, n + 1)
    V = zeros(Int, m + 1, n + 1)

    H[1, 1] = V[1, 1] = INF

    for i in 2:m+1
        M[i, 1] = INF
        H[i, 1] = INF
        V[i, 1] = gap_open + (i - 2) * gap_extend
    end
    for j in 2:n+1
        M[1, j] = INF
        H[1, j] = gap_open + (j - 2) * gap_extend
        V[1, j] = INF
    end

    @inbounds for i in 2:m+1
        rc = ref[i-1]
        for j in 2:n+1
            dscore = rc == query[j-1] ? match_score : mismatch
            M[i, j] = dscore + max(M[i-1, j-1], V[i-1, j-1], H[i-1, j-1])
            V[i, j] = max(H[i-1, j] + gap_open, V[i-1, j] + gap_extend, M[i-1, j] + gap_open)
            H[i, j] = max(H[i, j-1] + gap_extend, V[i, j-1] + gap_open, M[i, j-1] + gap_open)
        end
    end

    optimal = max(M[m+1, n+1], V[m+1, n+1], H[m+1, n+1])
    state = if optimal == M[m+1, n+1]; :M
    elseif optimal == V[m+1, n+1]; :V
    else; :H
    end

    ref_row = Char[]
    query_row = Char[]
    i, j = m + 1, n + 1
    errors = 0

    while i > 1 || j > 1
        if state == :M
            push!(ref_row, ref[i-1])
            push!(query_row, query[j-1])
            ref[i-1] != query[j-1] && (errors += 1)
            ms, vs, hs = M[i-1, j-1], V[i-1, j-1], H[i-1, j-1]
            i -= 1; j -= 1
        elseif state == :V
            push!(ref_row, ref[i-1])
            push!(query_row, '-')
            errors += 1
            ms = M[i-1, j] + gap_open
            vs = V[i-1, j] + gap_extend
            hs = H[i-1, j] + gap_open
            i -= 1
        else
            push!(ref_row, '-')
            push!(query_row, query[j-1])
            errors += 1
            ms = M[i, j-1] + gap_open
            vs = V[i, j-1] + gap_open
            hs = H[i, j-1] + gap_extend
            j -= 1
        end
        state = hs > ms && hs > vs ? :H : (vs > ms ? :V : :M)
    end

    Alignment(String(reverse(ref_row)), String(reverse(query_row)), optimal, errors)
end

"""
    describe_nt_change(ref::String, query::String) -> String

Describe nucleotide changes between ref and query using HGVS-like notation.
"""
function describe_nt_change(ref::AbstractString, query::AbstractString)
    aln = align_affine(ref, query)
    changes = String[]
    idx = 1
    i = 1
    while i <= length(aln.ref_row)
        c1, c2 = aln.ref_row[i], aln.query_row[i]
        if c1 == c2
            idx += 1
            i += 1
        elseif c1 == '-'  # insertion
            inserted = IOBuffer()
            while i <= length(aln.ref_row) && aln.ref_row[i] == '-'
                write(inserted, aln.query_row[i])
                i += 1
            end
            push!(changes, "$(idx-1)_$(idx)ins$(String(take!(inserted)))")
        elseif c2 == '-'  # deletion
            deleted = IOBuffer()
            start_idx = idx
            while i <= length(aln.ref_row) && aln.query_row[i] == '-'
                write(deleted, aln.ref_row[i])
                i += 1
                idx += 1
            end
            push!(changes, "$(start_idx)_$(idx-1)del$(String(take!(deleted)))")
        else  # substitution
            push!(changes, "$(idx)$(c1)>$(c2)")
            idx += 1
            i += 1
        end
    end
    join(changes, "; ")
end
