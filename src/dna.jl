# DNA sequence utilities — delegates to BioSequences for translation and complement

"""
    translate(seq) -> String

Translate nucleotide sequence to amino acids using the standard genetic code.
Ambiguous codons that resolve to a single amino acid are allowed; unresolvable ones become 'X'.
"""
function translate(seq::AbstractString)
    n = ncodeunits(seq) ÷ 3
    n == 0 && return ""
    dna = BioSequences.LongDNA{4}(uppercase(SubString(seq, 1, 3n)))
    String(BioSequences.translate(dna; allow_ambiguous_codons=true))
end

"""
    reverse_complement(seq) -> String

Return the reverse complement of a DNA sequence.
"""
function reverse_complement(seq::AbstractString)
    String(BioSequences.reverse_complement(BioSequences.LongDNA{4}(seq)))
end

"""
    has_stop(seq) -> Bool

Check whether sequence contains an internal stop codon. Incomplete trailing codon allowed.
"""
function has_stop(seq::AbstractString)
    n = ncodeunits(seq) ÷ 3
    n == 0 && return false
    occursin('*', translate(SubString(seq, 1, 3n)))
end

# ─── Thread-local edit distance buffers ───

struct EditDistanceBuffers
    prev::Vector{Int}
    curr::Vector{Int}
end

EditDistanceBuffers() = EditDistanceBuffers(Int[], Int[])

function ensure_capacity!(buf::EditDistanceBuffers, n::Int)
    needed = n + 1
    length(buf.prev) < needed && resize!(buf.prev, needed)
    length(buf.curr) < needed && resize!(buf.curr, needed)
    buf
end

const EDIT_BUFFERS = EditDistanceBuffers[]
const EDIT_BUFFERS_LOCK = ReentrantLock()

function get_edit_buffers()
    tid = Threads.threadid()
    if tid <= length(EDIT_BUFFERS)
        return @inbounds EDIT_BUFFERS[tid]
    end
    lock(EDIT_BUFFERS_LOCK) do
        while length(EDIT_BUFFERS) < tid
            push!(EDIT_BUFFERS, EditDistanceBuffers())
        end
    end
    @inbounds EDIT_BUFFERS[tid]
end

"""
    edit_distance(s, t; maxdiff) -> Int

Levenshtein distance with optional early termination at maxdiff+1.
Uses thread-local pre-allocated buffers for zero-allocation hot loops.
"""
function edit_distance(s::AbstractString, t::AbstractString; maxdiff::Int=typemax(Int))
    sv = codeunits(s)
    tv = codeunits(t)
    m, n = length(sv), length(tv)
    unbounded = (maxdiff == typemax(Int))
    if !unbounded && abs(m - n) > maxdiff
        return maxdiff + 1
    end

    buf = ensure_capacity!(get_edit_buffers(), n)
    prev = buf.prev
    curr = buf.curr

    @inbounds for j in 0:n
        prev[j+1] = j
    end

    @inbounds for i in 1:m
        curr[1] = i
        row_min = i
        si = sv[i]
        for j in 1:n
            cost = si == tv[j] ? 0 : 1
            curr[j+1] = min(curr[j] + 1, prev[j+1] + 1, prev[j] + cost)
            row_min = min(row_min, curr[j+1])
        end
        if !unbounded && row_min > maxdiff
            return maxdiff + 1
        end
        prev, curr = curr, prev
    end
    unbounded ? prev[n+1] : min(prev[n+1], maxdiff + 1)
end

"""
    hamming_distance(s, t) -> Int

Hamming distance between two equal-length strings.
"""
function hamming_distance(s::AbstractString, t::AbstractString)
    sv, tv = codeunits(s), codeunits(t)
    length(sv) == length(tv) || error("Strings must have equal length for Hamming distance")
    d = 0
    @inbounds for i in eachindex(sv)
        d += sv[i] != tv[i]
    end
    d
end

"""
    sequence_hash(seq; digits=4) -> String

Fingerprint like "S1234" derived from sequence content.
"""
function sequence_hash(seq::AbstractString; digits::Int=4)
    h = bytes2hex(sha256(Vector{UInt8}(codeunits(seq))))
    n = parse(UInt64, h[end-3:end]; base=16)
    "S" * lpad(string(n % 10^digits), digits, '0')
end

"""
    unique_name(name, seq) -> String

Gene name with sequence-derived suffix. Strips existing `_S....` suffix first.
"""
function unique_name(name::AbstractString, seq::AbstractString)
    base = first(split(name, "_S"; limit=2))
    base * "_" * sequence_hash(seq)
end
