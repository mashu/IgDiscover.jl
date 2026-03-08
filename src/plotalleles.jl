# Allele usage analysis — compute co-expression matrices for V/D/J alleles
#
# Produces a cross-tabulation of gene usage counts (e.g. V gene × J allele),
# exported as TSV for external plotting or visualization.

"""
    AlleleFilterCriteria(; d_evalue=1e-4, d_coverage=65.0)

Filter criteria for D gene quality when D is on either axis.
"""
struct AlleleFilterCriteria
    d_evalue::Float64
    d_coverage::Float64
end

AlleleFilterCriteria(; d_evalue::Float64=1e-4, d_coverage::Float64=65.0) =
    AlleleFilterCriteria(d_evalue, d_coverage)

"""
    AlleleUsage(; x_gene=:V, y_gene=:J, filter_criteria=AlleleFilterCriteria())

Callable struct that computes allele co-expression matrices.

```julia
au = AlleleUsage(x_gene=:V, y_gene=:J)
matrix = au(table, ["IGHJ4*02", "IGHJ6*02"])
```
"""
struct AlleleUsage
    x_gene::Symbol
    y_gene::Symbol
    filter_criteria::AlleleFilterCriteria
end

AlleleUsage(; x_gene::Symbol=:V, y_gene::Symbol=:J,
              filter_criteria::AlleleFilterCriteria=AlleleFilterCriteria()) =
    AlleleUsage(x_gene, y_gene, filter_criteria)

function gene_call_column(gene::Symbol)
    gene === :V && return :v_call
    gene === :D && return :d_call
    gene === :J && return :j_call
    error("Unknown gene type: $gene (expected :V, :D, or :J)")
end

"""
    apply_filters(table, au::AlleleUsage) -> DataFrame

Apply V/J error and D quality filters based on which genes are on the axes.
"""
function apply_filters(table::DataFrame, au::AlleleUsage)
    df = copy(table)
    if (au.x_gene === :V || au.y_gene === :V) && hasproperty(df, :V_errors)
        df = df[df.V_errors .== 0, :]
    end
    if (au.x_gene === :J || au.y_gene === :J) && hasproperty(df, :J_errors)
        df = df[df.J_errors .== 0, :]
    end
    if au.x_gene === :D || au.y_gene === :D
        fc = au.filter_criteria
        hasproperty(df, :d_support) && (df = df[coalesce.(df.d_support, Inf) .<= fc.d_evalue, :])
        hasproperty(df, :D_covered) && (df = df[coalesce.(df.D_covered, 0.0) .>= fc.d_coverage, :])
        hasproperty(df, :D_errors) && (df = df[df.D_errors .== 0, :])
    end
    df
end

"""
    (au::AlleleUsage)(table, alleles; database_names=nothing, order_names=nothing) -> DataFrame

Compute the allele co-expression matrix. Rows are x_gene assignments,
columns are the specified alleles, cells are co-occurrence counts.
"""
function (au::AlleleUsage)(table::DataFrame, alleles::Vector{String};
                           database_names::Union{Nothing,Vector{String}}=nothing,
                           order_names::Union{Nothing,Vector{String}}=nothing)
    df = apply_filters(table, au)
    x_col = gene_call_column(au.x_gene)
    y_col = gene_call_column(au.y_gene)

    x_vals = String.(coalesce.(df[!, x_col], ""))
    y_vals = String.(coalesce.(df[!, y_col], ""))

    counts = Dict{Tuple{String,String}, Int}()
    for i in eachindex(x_vals)
        (isempty(x_vals[i]) || isempty(y_vals[i])) && continue
        key = (x_vals[i], y_vals[i])
        counts[key] = get(counts, key, 0) + 1
    end

    x_names = sort(unique(filter(!isempty, x_vals)))
    database_names !== nothing && (x_names = filter(n -> n in Set(x_names), database_names))

    if order_names !== nothing
        gene_order = Dict(first(split(name, '*')) => idx for (idx, name) in enumerate(order_names))
        sort!(x_names; by=full_name -> begin
            base = first(split(full_name, '*'))
            allele_num = let parts = split(full_name, '*')
                length(parts) >= 2 ? something(tryparse(Int, parts[2]), 0) : 0
            end
            (get(gene_order, base, 1_000_000), allele_num)
        end)
    end

    all_y = Set(y_vals)
    for a in alleles
        a in all_y || @warn "Allele $a not expressed in this dataset"
    end

    result = DataFrame(gene = x_names)
    for a in alleles
        result[!, Symbol(a)] = [get(counts, (x, a), 0) for x in x_names]
    end
    result
end

"""
    allele_expression_counts(table, au::AlleleUsage) -> DataFrame

Per-gene expression counts summed across the y_gene axis.
"""
function allele_expression_counts(table::DataFrame, au::AlleleUsage)
    df = apply_filters(table, au)
    y_col = gene_call_column(au.y_gene)
    y_vals = String.(coalesce.(df[!, y_col], ""))
    counts = tallies(filter(!isempty, y_vals))
    genes = sort(collect(keys(counts)))
    DataFrame(gene=genes, count=[counts[g] for g in genes])
end
