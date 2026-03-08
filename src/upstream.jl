# Upstream sequence analysis — compute UTR and leader consensus per gene

"""
    UpstreamParams(;
        max_v_errors=1.0, max_fr1_errors=NaN, max_cdr1_errors=NaN,
        min_consensus_size=1, consensus_threshold=75.0,
        discard_ambiguous=false, part=:UTR_leader)

Parameters for upstream analysis. `part` selects :UTR, :leader, or :UTR_leader.
NaN means "do not filter on this field".
"""
struct UpstreamParams
    max_v_errors::Float64
    max_fr1_errors::Float64
    max_cdr1_errors::Float64
    min_consensus_size::Int
    consensus_threshold::Float64
    discard_ambiguous::Bool
    part::Symbol
end

function UpstreamParams(;
    max_v_errors::Float64=1.0,
    max_fr1_errors::Float64=NaN,
    max_cdr1_errors::Float64=NaN,
    min_consensus_size::Int=1,
    consensus_threshold::Float64=75.0,
    discard_ambiguous::Bool=false,
    part::Symbol=:UTR_leader,
)
    part in (:UTR, :leader, :UTR_leader) ||
        error("part must be :UTR, :leader, or :UTR_leader; got $part")
    UpstreamParams(max_v_errors, max_fr1_errors, max_cdr1_errors,
                   min_consensus_size, consensus_threshold, discard_ambiguous, part)
end

"""
    apply_filters(table, p::UpstreamParams) -> DataFrame

Apply SHM-based filters for upstream analysis.
"""
function apply_filters(table::DataFrame, p::UpstreamParams)
    df = copy(table)
    !isnan(p.max_v_errors) && hasproperty(df, :V_SHM) &&
        (df = df[coalesce.(df.V_SHM, Inf) .<= p.max_v_errors, :])
    !isnan(p.max_fr1_errors) && hasproperty(df, :FR1_SHM) &&
        (df = df[coalesce.(df.FR1_SHM, Inf) .<= p.max_fr1_errors, :])
    !isnan(p.max_cdr1_errors) && hasproperty(df, :CDR1_SHM) &&
        (df = df[coalesce.(df.CDR1_SHM, Inf) .<= p.max_cdr1_errors, :])
    df
end

"""
    upstream_target_column(df, part) -> Vector{String}

Build the target sequence column for the specified `part`.
"""
function upstream_target_column(df::DataFrame, part::Symbol)
    ensure_column!(df, :UTR, "")
    ensure_column!(df, :leader, "")
    utrs = String.(df.UTR)
    leaders = String.(df.leader)
    part === :UTR ? utrs : part === :leader ? leaders : utrs .* leaders
end

"""
    UpstreamAnalyzer(params::UpstreamParams)

Callable struct computing per-gene upstream consensus sequences.

```julia
analyzer = UpstreamAnalyzer(UpstreamParams(max_v_errors=1.0))
records = analyzer(table)  # returns Vector{FastaRecord}
```
"""
struct UpstreamAnalyzer
    params::UpstreamParams
end

UpstreamAnalyzer(; kwargs...) = UpstreamAnalyzer(UpstreamParams(; kwargs...))

"""
    (analyzer::UpstreamAnalyzer)(table) -> Vector{FastaRecord}
"""
function (analyzer::UpstreamAnalyzer)(table::DataFrame)
    p = analyzer.params
    df = apply_filters(table, p)
    target_col = upstream_target_column(df, p.part)
    keep = .!isempty.(target_col)
    df = df[keep, :]
    target_col = target_col[keep]

    ensure_column!(df, :v_call, "")
    v_calls = String.(df.v_call)

    records = FastaRecord[]
    n_ambiguous = 0

    for gene in sort(unique(filter(!isempty, v_calls)))
        mask = v_calls .== gene
        gene_targets = target_col[mask]
        gene_utrs = hasproperty(df, :UTR) ? String.(coalesce.(df[mask, :UTR], "")) : gene_targets

        length(gene_targets) < p.min_consensus_size && continue

        sequences = if p.part === :leader
            gene_targets
        else
            utr_lengths = length.(gene_utrs)
            sorted_lengths = sort(utr_lengths; rev=true)
            length_threshold = sorted_lengths[min(10, length(sorted_lengths))]
            gene_targets[utr_lengths .>= length_threshold]
        end

        isempty(sequences) && continue

        cons = if length(sequences) == 1
            sequences[1]
        else
            lstrip(iterative_consensus(sequences;
                program="muscle-medium",
                threshold=p.consensus_threshold / 100), 'N')
        end

        n_count = count(==('N'), cons)
        n_count > 0 && (n_ambiguous += 1)
        (n_count > 0 && p.discard_ambiguous) && continue

        push!(records, FastaRecord(gene, cons))
    end

    @info "Wrote consensus for $(length(records)) genes ($(n_ambiguous) with ambiguous bases)"
    records
end
