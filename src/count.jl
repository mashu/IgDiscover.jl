# Expression count computation — count per-gene assignment frequencies

"""
    ExpressionCounter(; gene=:V, d_evalue=nothing, d_coverage=nothing,
                        d_errors=nothing, allele_ratio=nothing)

Callable struct computing gene expression counts. NaN sentinel = "not set".

```julia
counter = ExpressionCounter(gene=:V)
counts = counter(table; gene_names=["IGHV1*01"])
```
"""
struct ExpressionCounter
    gene::Symbol
    d_evalue::Float64
    d_coverage::Float64
    d_errors::Float64
    allele_ratio::Float64
end

function ExpressionCounter(;
    gene::Symbol=:V,
    d_evalue::Union{Nothing,Float64}=nothing,
    d_coverage::Union{Nothing,Float64}=nothing,
    d_errors::Union{Nothing,Int}=nothing,
    allele_ratio::Union{Nothing,Float64}=nothing,
)
    gene in (:V, :D, :J) || error("gene must be :V, :D, or :J; got $gene")
    de = d_evalue !== nothing ? d_evalue : (gene === :D ? 1e-4 : NaN)
    dc = d_coverage !== nothing ? d_coverage : (gene === :D ? 70.0 : NaN)
    dd = d_errors !== nothing ? Float64(d_errors) : NaN
    ar = allele_ratio !== nothing ? allele_ratio : NaN
    ExpressionCounter(gene, de, dc, dd, ar)
end

"""
    apply_filters(table, counter::ExpressionCounter) -> DataFrame

Apply D gene quality filters for expression counting.
"""
function apply_filters(table::DataFrame, counter::ExpressionCounter)
    df = table
    !isnan(counter.d_evalue) && hasproperty(df, :d_support) &&
        (df = df[coalesce.(df.d_support, Inf) .<= counter.d_evalue, :])
    !isnan(counter.d_coverage) && hasproperty(df, :D_covered) &&
        (df = df[coalesce.(df.D_covered, 0.0) .>= counter.d_coverage, :])
    !isnan(counter.d_errors) && hasproperty(df, :D_errors) &&
        (df = df[coalesce.(df.D_errors, typemax(Int)) .<= Int(counter.d_errors), :])
    df
end

"""
    (counter::ExpressionCounter)(table; gene_names=nothing) -> DataFrame

Count gene expression. Returns DataFrame with `gene` and `count` columns.
"""
function (counter::ExpressionCounter)(table::DataFrame;
            gene_names::Union{Nothing,Vector{String}}=nothing)
    df = apply_filters(table, counter)
    col = gene_call_column(counter.gene)
    hasproperty(df, col) || error("Table missing column $col")

    gene_vals = String.(coalesce.(df[!, col], ""))
    counts = tallies(filter(!isempty, gene_vals))

    if gene_names !== nothing
        genes = gene_names
        count_vals = [get(counts, g, 0) for g in genes]
    else
        genes = sort(collect(keys(counts)); by=natural_sort_key)
        count_vals = [counts[g] for g in genes]
    end

    result = DataFrame(gene=genes, count=count_vals)
    !isnan(counter.allele_ratio) && return filter_by_allele_ratio(result, counter.allele_ratio)
    result
end

"""
    filter_by_allele_ratio(counts, ratio) -> DataFrame

Remove alleles below `ratio` × best allele count within the same gene family.
"""
function filter_by_allele_ratio(counts::DataFrame, ratio::Float64)
    families = gene_family.(counts.gene)
    best = Dict{String,Int}()
    for (fam, cnt) in zip(families, counts.count)
        best[fam] = max(get(best, fam, 0), cnt)
    end
    mask = [let fb = get(best, families[i], 0)
                fb == 0 || counts.count[i] / fb >= ratio
            end for i in 1:nrow(counts)]
    counts[mask, :]
end

"""
    natural_sort_key(s) -> Vector

Sorting key for natural order: "file2" < "file10".
"""
function natural_sort_key(s::AbstractString)
    parts = split(s, r"([0-9]+)"; keepempty=false)
    map(parts) do p
        n = tryparse(Int, p)
        n !== nothing ? (0, n, "") : (1, 0, lowercase(p))
    end
end
