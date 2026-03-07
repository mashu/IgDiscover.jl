# Shared utilities — pure functions used across modules

"""
    tallies(v) -> Dict{eltype(v), Int}

Count occurrences of each element in `v`.
"""
function tallies(v::AbstractVector)
    d = Dict{eltype(v),Int}()
    for x in v
        d[x] = get(d, x, 0) + 1
    end
    d
end

"""
    most_common(v) -> (element, count)

Return the most frequent element and its count.
"""
function most_common(v::AbstractVector)
    t = tallies(v)
    isempty(t) && return (first(v), 0)
    first(sort!(collect(t); by=last, rev=true))
end

"""
    safe_divide(a, b) -> Float64

Division that returns 0.0 when `b == 0`.
"""
safe_divide(a, b) = b == 0 ? 0.0 : Float64(a) / Float64(b)

"""
    is_same_gene(name1, name2) -> Bool

Check if two allele names belong to the same gene (same prefix before '*').
"""
function is_same_gene(name1::AbstractString, name2::AbstractString)
    g1 = occursin('*', name1) ? first(split(name1, '*')) : name1
    g2 = occursin('*', name2) ? first(split(name2, '*')) : name2
    g1 == g2
end

# ─── Generic TOML → struct constructor ───

"""
    from_toml(T, d::Dict) -> T

Construct a struct `T` from a TOML dictionary by matching field names to keys.
Each field is converted to its declared type. Nested structs are constructed recursively
if the corresponding value is a `Dict`.
"""
function from_toml(::Type{T}, d::Dict) where T
    args = ntuple(fieldcount(T)) do i
        fname = String(fieldname(T, i))
        ftype = fieldtype(T, i)
        raw = d[fname]
        if raw isa Dict && !(ftype <: Dict)
            from_toml(ftype, raw)
        else
            convert(ftype, raw)
        end
    end
    T(args...)
end
