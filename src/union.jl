# Union — compute union of sequences across multiple FASTA files
#
# Merges sequences where one is a prefix of the other, keeping the longer one.
# This is useful for combining gene databases from different sources.

"""
    SequenceUnion()

Callable struct that merges FASTA records from multiple files, combining
sequences where one is a prefix of the other into a single entry (keeping
the longer sequence and the first name seen).

```julia
union = SequenceUnion()
merged = union(["db1.fasta", "db2.fasta"])
write_fasta("combined.fasta", merged)
```
"""
struct SequenceUnion end

struct NamedSequence
    name::String
    sequence::String
end

"""
    (::SequenceUnion)(paths) -> Vector{FastaRecord}

Read all FASTA files and merge sequences that share a prefix relationship.
"""
function (::SequenceUnion)(paths::Vector{String})
    items = NamedSequence[]
    n_read = 0
    for path in paths
        for rec in read_fasta(path)
            items = merge_into(items, NamedSequence(rec.name, uppercase(rec.sequence)))
            n_read += 1
        end
        @info "Read sequences from $path"
    end
    @info "Read $n_read sequences from $(length(paths)) files, merged to $(length(items))"
    sorted = sort(items; by=ns -> ns.name)
    [FastaRecord(ns.name, ns.sequence) for ns in sorted]
end

"""
    merge_into(items, new_item) -> Vector{NamedSequence}

Attempt to merge `new_item` with existing items. Two sequences merge when one
is a prefix of the other; the longer one is kept with the first name.
"""
function merge_into(items::Vector{NamedSequence}, item::NamedSequence)
    result = NamedSequence[]
    merged = item
    for existing in items
        m = try_merge(existing, merged)
        if m === nothing
            push!(result, existing)
        else
            merged = m
        end
    end
    push!(result, merged)
    result
end

"""
    try_merge(a, b) -> Union{Nothing, NamedSequence}

Merge two named sequences if one is a prefix of the other.
Returns the longer one (with `a`'s name if lengths match), or `nothing` if unrelated.
"""
function try_merge(a::NamedSequence, b::NamedSequence)
    sa = a.sequence
    sb = b.sequence
    # Pad shorter to same length to ignore end-gap differences
    la, lb = length(sa), length(sb)
    if la < lb
        sa = sa * sb[la+1:lb]
    elseif lb < la
        sb = sb * sa[lb+1:la]
    end
    sa == sb || return nothing
    length(a.sequence) >= length(b.sequence) ? a : NamedSequence(a.name, b.sequence)
end
