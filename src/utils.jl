# Shared utilities — pure functions used across modules

"""
    tallies(v) -> Dict{T, Int}

Count occurrences of each element in `v`. Uses parametric type `T` from `AbstractVector{T}`.
"""
function tallies(v::AbstractVector{T}) where T
    d = Dict{T,Int}()
    for x in v
        d[x] = get(d, x, 0) + 1
    end
    d
end

"""
    most_common(v) -> (element, count)

Return the most frequent element and its count.
"""
function most_common(v::AbstractVector{T}) where T
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
    gene_family(name) -> String

Extract gene family (everything before '*') from an allele name.
Returns the full name if no '*' is present.
"""
gene_family(name::AbstractString) =
    occursin('*', name) ? first(split(name, '*')) : String(name)

"""
    is_same_gene(name1, name2) -> Bool

Check if two allele names belong to the same gene (same prefix before '*').
"""
is_same_gene(name1::AbstractString, name2::AbstractString) =
    gene_family(name1) == gene_family(name2)

# ─── DataFrame column initialization ───
#
# These eliminate the repeated hasproperty/fill/coalesce boilerplate
# that appears in augment.jl, discovery.jl, jdiscovery.jl, and clonotypes.jl.

"""
    ensure_column!(df, col, default::String)
    ensure_column!(df, col, default::Int)
    ensure_column!(df, col, default::Float64)

Ensure `col` exists in `df` with the correct element type, filling missing values with `default`.
Dispatches on the type of `default`.
"""
function ensure_column!(df::DataFrame, col::Symbol, default::String)
    n = nrow(df)
    hasproperty(df, col) || (df[!, col] = fill(default, n))
    df[!, col] = Vector{String}(coalesce.(df[!, col], default))
end

function ensure_column!(df::DataFrame, col::Symbol, default::Int)
    n = nrow(df)
    hasproperty(df, col) || (df[!, col] = fill(default, n))
    df[!, col] = Vector{Int}(coalesce.(df[!, col], default))
end

function ensure_column!(df::DataFrame, col::Symbol, default::Float64)
    n = nrow(df)
    hasproperty(df, col) || (df[!, col] = fill(default, n))
    df[!, col] = Vector{Float64}(coalesce.(df[!, col], default))
end

# ─── Generic TOML → struct constructor ───
#
# Dispatch on value type: nested Dict → recurse; scalar → convert.
# No runtime type checks; compile-time dispatch only.

from_toml_field(::Type{T}, raw::Dict{K,V}) where {T,K,V} = from_toml(T, raw)
from_toml_field(::Type{T}, raw) where T = convert(T, raw)

"""
    from_toml(T, d::Dict) -> T

Construct a struct `T` from a TOML dictionary by matching field names to keys.
Each field is converted to its declared type. Nested structs are constructed recursively
via dispatch when the value is a `Dict`.
"""
function from_toml(::Type{T}, d::Dict{K,V}) where {T,K,V}
    args = ntuple(fieldcount(T)) do i
        fname = String(fieldname(T, i))
        ftype = fieldtype(T, i)
        raw = d[fname]
        from_toml_field(ftype, raw)
    end
    T(args...)
end
