# Database comparison — compare two FASTA files based on sequences

"""
    DatabaseDiff

Categorized result of comparing two FASTA databases.
"""
struct DatabaseDiff
    identical::Vector{Tuple{FastaRecord, FastaRecord}}
    similar::Vector{Tuple{FastaRecord, FastaRecord}}
    only_a::Vector{FastaRecord}
    only_b::Vector{FastaRecord}
    renamed::Vector{Tuple{FastaRecord, FastaRecord}}
end

"""
    exit_code(diff::DatabaseDiff) -> Int

0 if identical (name differences allowed), 1 if lost/gained/mutated.
"""
exit_code(diff::DatabaseDiff) =
    (isempty(diff.only_a) && isempty(diff.only_b) && isempty(diff.similar)) ? 0 : 1

"""
    DatabaseComparator(max_cost=20)

Callable struct comparing two FASTA record sets. Pairs identical sequences
first, then uses greedy minimum-cost matching for remaining records.

```julia
cmp = DatabaseComparator()
diff = cmp(read_fasta("a.fasta"), read_fasta("b.fasta"))
```
"""
struct DatabaseComparator
    max_cost::Int
end
DatabaseComparator() = DatabaseComparator(20)

function (cmp::DatabaseComparator)(a_records::Vector{FastaRecord}, b_records::Vector{FastaRecord})
    identical, a_rest, b_rest = pair_identical(a_records, b_records)
    similar_pairs, a_only, b_only = pair_similar(a_rest, b_rest, cmp.max_cost)
    renamed = filter(((a, b),) -> a.name != b.name, identical)
    DatabaseDiff(identical, similar_pairs, a_only, b_only, renamed)
end

# ─── Comparison cost ───

function comparison_cost(a::FastaRecord, b::FastaRecord)
    la, lb = length(a.sequence), length(b.sequence)
    l = min(la, lb)
    length_diff = max(la, lb) - l
    dp = hamming_distance(a.sequence[1:l], b.sequence[1:l])
    ds = hamming_distance(a.sequence[end-l+1:end], b.sequence[end-l+1:end])
    5 * min(dp, ds) + length_diff
end

# ─── Phase 1: pair identical ───

function pair_identical(a_records::Vector{FastaRecord}, b_records::Vector{FastaRecord})
    b_map = Dict(b.sequence => b for b in b_records)
    identical = Tuple{FastaRecord, FastaRecord}[]
    a_rest = FastaRecord[]
    for a in a_records
        if haskey(b_map, a.sequence)
            push!(identical, (a, b_map[a.sequence]))
            delete!(b_map, a.sequence)
        else
            push!(a_rest, a)
        end
    end
    (identical, a_rest, collect(values(b_map)))
end

# ─── Phase 2: greedy pairing ───

function pair_similar(a_rest::Vector{FastaRecord}, b_rest::Vector{FastaRecord}, max_cost::Int)
    m, n = length(a_rest), length(b_rest)
    (m == 0 || n == 0) && return (Tuple{FastaRecord,FastaRecord}[], a_rest, b_rest)

    assignments = [(comparison_cost(a_rest[i], b_rest[j]), i, j)
                   for i in 1:m for j in 1:n]
    sort!(assignments; by=first)

    a_matched = falses(m)
    b_matched = falses(n)
    similar_pairs = Tuple{FastaRecord, FastaRecord}[]

    for (c, i, j) in assignments
        c > max_cost && break
        (a_matched[i] || b_matched[j]) && continue
        push!(similar_pairs, (a_rest[i], b_rest[j]))
        a_matched[i] = true
        b_matched[j] = true
    end
    (similar_pairs, a_rest[.!a_matched], b_rest[.!b_matched])
end

# ─── Formatting ───

"""
    format_diff(io, diff; color=false)

Write a human-readable comparison report.
"""
function format_diff(io::IO, diff::DatabaseDiff; color::Bool=false)
    n_a = length(diff.identical) + length(diff.only_a) + length(diff.similar)
    n_b = length(diff.identical) + length(diff.only_b) + length(diff.similar)
    n_ident = length(diff.identical) - length(diff.renamed)
    println(io, "$n_a vs $n_b records. $(length(diff.only_a)) lost, ",
            "$(length(diff.only_b)) gained, $n_ident identical, ",
            "$(length(diff.renamed)) different name, $(length(diff.similar)) similar")

    if !isempty(diff.only_a)
        println(io, "\n## Only in A")
        for r in diff.only_a; println(io, "- ", r.name); end
    end
    if !isempty(diff.only_b)
        println(io, "\n## Only in B")
        for r in diff.only_b; println(io, "+ ", r.name); end
    end
    if !isempty(diff.renamed)
        println(io, "\n## Different name (sequence identical)")
        for (a, b) in diff.renamed; println(io, "= ", a.name, " -- ", b.name); end
    end
    if !isempty(diff.similar)
        println(io, "\n## Similar")
        for (a, b) in diff.similar
            format_similar_pair(io, a, b)
        end
    end
end

function format_similar_pair(io::IO, a::FastaRecord, b::FastaRecord)
    println(io, "~ ", a.name, " -- ", b.name)
    la, lb = length(a.sequence), length(b.sequence)
    l = min(la, lb)
    dp = hamming_distance(a.sequence[1:l], b.sequence[1:l])
    ds = hamming_distance(a.sequence[end-l+1:end], b.sequence[end-l+1:end])
    ac, bc = dp <= ds ? (a.sequence[1:l], b.sequence[1:l]) :
                        (a.sequence[end-l+1:end], b.sequence[end-l+1:end])
    edits = [ac_c == bc_c ? string(ac_c) : "{$ac_c → $bc_c}"
             for (ac_c, bc_c) in zip(ac, bc)]
    println(io, join(edits))
    println(io)
end

# ─── Duplicate detection ───

"""
    check_duplicate_names(records) -> Vector{String}
"""
function check_duplicate_names(records::Vector{FastaRecord})
    seen = Set{String}()
    dups = String[]
    for r in records
        r.name in seen && push!(dups, r.name)
        push!(seen, r.name)
    end
    dups
end

"""
    check_duplicate_sequences(records) -> Vector{Tuple{String,String}}
"""
function check_duplicate_sequences(records::Vector{FastaRecord})
    seen = Dict{String, String}()
    dups = Tuple{String,String}[]
    for r in records
        haskey(seen, r.sequence) ? push!(dups, (r.name, seen[r.sequence])) :
                                   (seen[r.sequence] = r.name)
    end
    dups
end
