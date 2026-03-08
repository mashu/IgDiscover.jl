# Haplotype analysis — determine haplotypes based on allele co-occurrences
#
# Uses heterozygous genes (two alleles with sufficient expression) to phase
# the remaining genes into two haplotype groups via co-expression counting.

const HETEROZYGOUS_THRESHOLD = 0.1
const EXPRESSED_RATIO = 0.1

# ─── Types ───

"""
    HaplotypeEntry(hap1, hap2, classification, counts)

A single gene's haplotype assignment. `hap1`/`hap2` are allele names
(empty string if absent on that haplotype). `classification` is one of
:homozygous, :heterozygous, :deletion, :duplication, :unknown.
"""
struct HaplotypeEntry
    hap1::String
    hap2::String
    classification::Symbol
    count1::Int
    count2::Int
end

"""
    HaplotypePair(entries, gene_type, het1, het2)

Haplotype pair for a single gene type (V/D/J), phased using the heterozygous
gene with alleles `het1` and `het2`.
"""
mutable struct HaplotypePair
    entries::Vector{HaplotypeEntry}
    gene_type::Char
    het1::String
    het2::String
end

function switch_haplotype!(hp::HaplotypePair)
    hp.het1, hp.het2 = hp.het2, hp.het1
    hp.entries = [HaplotypeEntry(e.hap2, e.hap1, e.classification, e.count2, e.count1)
                  for e in hp.entries]
end

function sort_haplotype!(hp::HaplotypePair, gene_order::Dict{String,Int})
    sort!(hp.entries; by=e -> begin
        name = isempty(e.hap1) ? e.hap2 : e.hap1
        base, _, allele_str = split_allele(name)
        allele_num = something(tryparse(Int, allele_str), 999)
        get(gene_order, base, 1_000_000) * 1000 + allele_num
    end)
end

function format_tsv(hp::HaplotypePair; header::Bool=true)
    lines = String[]
    header && push!(lines, join(["haplotype1", "haplotype2", "type", "count1", "count2"], '\t'))
    push!(lines, "# $(hp.gene_type) haplotype from $(hp.het1) and $(hp.het2)")
    for e in hp.entries
        push!(lines, join([e.hap1, e.hap2, string(e.classification),
                           string(e.count1), string(e.count2)], '\t'))
    end
    join(lines, '\n') * '\n'
end

# ─── Helpers ───

function split_allele(name::AbstractString)
    idx = findfirst('*', name)
    idx === nothing && return (name, "", "")
    (name[1:idx-1], "*", name[idx+1:end])
end

# ─── Expression counting ───

struct AlleleExpression
    gene::String
    allele::String
    name::String
    count::Int
end

"""
    expression_by_gene(table, gene_type) -> Vector{Vector{AlleleExpression}}

Group allele expressions by gene, filtering low-expression alleles.
Each inner vector contains alleles for one gene.
"""
function expression_by_gene(table::DataFrame, gene_type::Char)
    col = Symbol(lowercase(string(gene_type)), "_call")
    hasproperty(table, col) || return Vector{AlleleExpression}[]

    gene_vals = String.(coalesce.(table[!, col], ""))
    counts = tallies(filter(!isempty, gene_vals))

    # Parse into gene/allele pairs
    entries = AlleleExpression[]
    for (name, cnt) in counts
        base, _, allele = split_allele(name)
        push!(entries, AlleleExpression(base, allele, name, cnt))
    end

    # Group by gene family
    groups = Dict{String, Vector{AlleleExpression}}()
    for e in entries
        push!(get!(groups, e.gene, AlleleExpression[]), e)
    end

    result = Vector{AlleleExpression}[]
    for (_, alleles) in sort!(collect(groups); by=first)
        max_expr = maximum(a.count for a in alleles)
        filtered = filter(a -> a.count >= HETEROZYGOUS_THRESHOLD * max_expr, alleles)
        sort!(filtered; by=a -> a.count, rev=true)
        push!(result, filtered)
    end
    result
end

# ─── Co-expression matrix ───

function compute_coexpressions(table::DataFrame, het_gene::Char, target_gene::Char)
    col1 = Symbol(lowercase(string(het_gene)), "_call")
    col2 = Symbol(lowercase(string(target_gene)), "_call")
    v1 = String.(coalesce.(table[!, col1], ""))
    v2 = String.(coalesce.(table[!, col2], ""))

    counts = Dict{Tuple{String,String}, Int}()
    for i in eachindex(v1)
        (isempty(v1[i]) || isempty(v2[i])) && continue
        key = (v1[i], v2[i])
        counts[key] = get(counts, key, 0) + 1
    end
    counts
end

# ─── Co-occurrence classification ───

function classify_cooccurrences(coexpressions::Dict{Tuple{String,String},Int},
                                het_alleles::Tuple{String,String},
                                target_groups::Vector{Vector{AlleleExpression}})
    entries = HaplotypeEntry[]

    for target_alleles in target_groups
        expressed_list = Vector{Tuple{Bool,Bool}}()
        allele_names = String[]
        allele_counts = Vector{Tuple{Int,Int}}()

        for ta in target_alleles
            c1 = get(coexpressions, (het_alleles[1], ta.name), 0)
            c2 = get(coexpressions, (het_alleles[2], ta.name), 0)
            total = c1 + c2 + 1  # avoid division by zero
            r1 = c1 / total >= EXPRESSED_RATIO
            r2 = c2 / total >= EXPRESSED_RATIO
            (r1 || r2) || continue
            push!(expressed_list, (r1, r2))
            push!(allele_names, ta.name)
            push!(allele_counts, (c1, c2))
        end

        if length(expressed_list) == 1
            r1, r2 = expressed_list[1]
            c1, c2 = allele_counts[1]
            n = allele_names[1]
            if r1 && !r2
                push!(entries, HaplotypeEntry(n, "", :deletion, c1, c2))
            elseif !r1 && r2
                push!(entries, HaplotypeEntry("", n, :deletion, c1, c2))
            else  # both true
                push!(entries, HaplotypeEntry(n, n, :homozygous, c1, c2))
            end
        elseif length(expressed_list) == 2 &&
               expressed_list == [(true, false), (false, true)]
            push!(entries, HaplotypeEntry(allele_names[1], allele_names[2],
                :heterozygous, allele_counts[1][1], allele_counts[2][2]))
        elseif length(expressed_list) == 2 &&
               expressed_list == [(false, true), (true, false)]
            push!(entries, HaplotypeEntry(allele_names[2], allele_names[1],
                :heterozygous, allele_counts[1][2], allele_counts[2][1]))
        else
            # Check for duplication pattern
            has_single = any(e -> e == (true, false) || e == (false, true), expressed_list)
            n_true = sum(count(identity, [e...]) for e in expressed_list)
            classification = (has_single && n_true > 2) ? :duplication : :unknown
            for (idx, (r1, r2)) in enumerate(expressed_list)
                push!(entries, HaplotypeEntry(
                    r1 ? allele_names[idx] : "",
                    r2 ? allele_names[idx] : "",
                    classification, allele_counts[idx]...))
            end
        end
    end
    entries
end

# ─── HaplotypeAnalyzer functor ───

"""
    HaplotypeAnalyzer(; d_evalue=1e-4, d_coverage=65.0)

Callable struct that determines haplotypes from a filtered assignment table.

```julia
analyzer = HaplotypeAnalyzer()
blocks = analyzer(table)
for block in blocks
    print(format_tsv(block; header=(block === first(blocks))))
end
```
"""
struct HaplotypeAnalyzer
    d_evalue::Float64
    d_coverage::Float64
end

HaplotypeAnalyzer(; d_evalue::Float64=1e-4, d_coverage::Float64=65.0) =
    HaplotypeAnalyzer(d_evalue, d_coverage)

"""
    apply_filters(table, ha::HaplotypeAnalyzer) -> DataFrame

Apply strict quality filters: V errors=0, J errors=0, D evalue/coverage.
"""
function apply_filters(table::DataFrame, ha::HaplotypeAnalyzer)
    df = copy(table)
    hasproperty(df, :V_errors) && (df = df[df.V_errors .== 0, :])
    hasproperty(df, :J_errors) && (df = df[df.J_errors .== 0, :])
    hasproperty(df, :d_support) && (df = df[coalesce.(df.d_support, Inf) .<= ha.d_evalue, :])
    hasproperty(df, :D_covered) && (df = df[coalesce.(df.D_covered, 0.0) .>= ha.d_coverage, :])
    hasproperty(df, :D_errors) && (df = df[df.D_errors .== 0, :])
    df
end

"""
    (ha::HaplotypeAnalyzer)(table; restrict_names=nothing, gene_order=nothing,
                            v_gene=nothing) -> Vector{HaplotypePair}
"""
function (ha::HaplotypeAnalyzer)(table::DataFrame;
            restrict_names::Union{Nothing,Set{String}}=nothing,
            gene_order::Union{Nothing,Vector{String}}=nothing,
            v_gene::Union{Nothing,String}=nothing)

    df = apply_filters(table, ha)
    if restrict_names !== nothing
        hasproperty(df, :v_call) &&
            (df = df[in.(String.(coalesce.(df.v_call, "")), Ref(restrict_names)), :])
    end

    # Compute expressions and find heterozygous genes
    expressions = Dict{Char, Vector{Vector{AlleleExpression}}}()
    het_expressions = Dict{Char, Vector{Union{Nothing,Vector{AlleleExpression}}}}()

    for gt in ('V', 'D', 'J')
        ex = expression_by_gene(df, gt)
        expressions[gt] = ex
        het_ex = [e for e in ex if length(e) == 2]
        if !isempty(het_ex)
            sort!(het_ex; by=g -> sum(a.count for a in g), rev=true)
            het_expressions[gt] = het_ex[1:min(5, length(het_ex))]
        else
            het_expressions[gt] = [nothing]
        end
    end

    # Override V gene if specified
    if v_gene !== nothing
        het_ex = [e for e in expressions['V'] if length(e) == 2]
        found = false
        for ex in het_ex
            if any(a -> a.name == v_gene || a.gene == v_gene, ex)
                het_expressions['V'] = [ex]
                found = true
                break
            end
        end
        found || @warn "Gene $v_gene not found among heterozygous V genes"
    end

    order_dict = gene_order !== nothing ?
        Dict(name => i for (i, name) in enumerate(gene_order)) : nothing

    # Try combinations of het J × het V
    best_blocks = HaplotypePair[]
    for het_j in het_expressions['J'], het_v in het_expressions['V']
        best_het = Dict{Char, Union{Nothing,Vector{AlleleExpression}}}(
            'V' => het_v,
            'D' => isempty(het_expressions['D']) ? nothing : het_expressions['D'][1],
            'J' => het_j,
        )

        blocks = HaplotypePair[]
        for (target_gt, het_gt) in [('J','V'), ('D','J'), ('V','J')]
            het_alleles = best_het[het_gt]
            het_alleles === nothing && continue
            coex = compute_coexpressions(df, het_gt, target_gt)
            het1, het2 = het_alleles[1].name, het_alleles[2].name
            entries = classify_cooccurrences(coex, (het1, het2), expressions[target_gt])
            block = HaplotypePair(entries, target_gt, het1, het2)
            order_dict !== nothing && sort_haplotype!(block, order_dict)
            push!(blocks, block)
        end

        # Check for duplication conflicts
        het_used = Set{String}()
        for (_, ha_val) in best_het
            ha_val === nothing && continue
            for a in ha_val; push!(het_used, a.name); end
        end

        has_conflict = false
        for block in blocks
            block.gene_type == 'D' && continue
            for e in block.entries
                for name in (e.hap1, e.hap2)
                    if name in het_used && e.classification != :heterozygous
                        has_conflict = true
                    end
                end
            end
        end

        if isempty(best_blocks) || !has_conflict
            best_blocks = blocks
            has_conflict || break
        end
    end

    # Phase across blocks: ensure J haplotype is consistent with V
    if length(best_blocks) == 3
        j_hap, d_hap, v_hap = best_blocks
        for e in j_hap.entries
            if (e.hap1, e.hap2) == (v_hap.het2, v_hap.het1)
                switch_haplotype!(j_hap)
                break
            end
        end
    end

    best_blocks
end
