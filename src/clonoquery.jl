# Clonoquery — query a reference table by clonotype from a query table
#
# For each query sequence, find all reference sequences that share the same
# clonotype (V gene, J gene, CDR3 length, CDR3 similarity within threshold).
# Optionally requires identical junction regions when a CDR3 core slice is given.

# ─── CDR3 similarity with optional junction constraint ───

"""
    CDR3CoreSlice(start, stop)

Define the non-junction region of CDR3 sequences. Positions before `start`
and after `stop` are junction regions; at least one junction must be identical
for two CDR3s to be considered similar.
"""
struct CDR3CoreSlice
    start::Int
    stop::Int
end

"""
    is_similar_with_junction(s, t, mismatches, core::Nothing) -> Bool

CDR3 similarity without junction constraint: same length and Hamming ≤ threshold.
When `mismatches` ∈ (0,1), it is a fraction of CDR3 length.
"""
function is_similar_with_junction(s::AbstractString, t::AbstractString,
                                  mismatches::Float64, ::Nothing)
    length(s) == length(t) || return false
    threshold = mismatches < 1.0 ? length(s) * mismatches : mismatches
    hamming_distance(s, t) <= threshold
end

"""
    is_similar_with_junction(s, t, mismatches, core::CDR3CoreSlice) -> Bool

CDR3 similarity with junction constraint: Hamming within threshold *and*
at least one junction region (before `core.start` or after `core.stop`)
is identical between `s` and `t`.
"""
function is_similar_with_junction(s::AbstractString, t::AbstractString,
                                  mismatches::Float64, core::CDR3CoreSlice)
    length(s) == length(t) || return false
    delta = core.start
    threshold = mismatches < 1.0 ? (length(s) - delta) * mismatches : mismatches
    hamming_distance(s, t) <= threshold || return false
    s[1:core.start] == t[1:core.start] || s[core.stop+1:end] == t[core.stop+1:end]
end

# ─── Query result ───

struct ClonoQueryResult
    query_id::String
    match_count::Int
    matches::DataFrame
end

# ─── Clonoquery functor ───

"""
    ClonoQuery(; mismatches=1.0, cdr3_core=nothing, use_aa=false, min_count=1)

Callable struct that queries a reference table by clonotype.

```julia
cq = ClonoQuery(mismatches=1.0)
results = cq(query_table, ref_table)
```
"""
struct ClonoQuery{C}
    mismatches::Float64
    cdr3_core::C
    use_aa::Bool
    min_count::Int
end

ClonoQuery(; mismatches::Float64=1.0, cdr3_core::Union{Nothing,CDR3CoreSlice}=nothing,
             use_aa::Bool=false, min_count::Int=1) =
    ClonoQuery(mismatches, cdr3_core, use_aa, min_count)

"""
    (cq::ClonoQuery)(query_table, ref_table) -> Vector{ClonoQueryResult}

For each query row, find all ref rows sharing the same clonotype.
"""
function (cq::ClonoQuery)(query_table::DataFrame, ref_table::DataFrame)
    cdr3_col = cq.use_aa ? :cdr3_aa : :cdr3
    for df in (query_table, ref_table)
        for col in (:v_call, :j_call, cdr3_col)
            hasproperty(df, col) || error("Table missing required column: $col")
        end
    end

    ref = cq.min_count > 1 && hasproperty(ref_table, :count) ?
        ref_table[ref_table.count .>= cq.min_count, :] : ref_table

    ref_cdr3 = String.(coalesce.(ref[!, cdr3_col], ""))
    ref_v    = String.(coalesce.(ref.v_call, ""))
    ref_j    = String.(coalesce.(ref.j_call, ""))
    ref_len  = length.(ref_cdr3)

    # Group reference by (v, j, cdr3_length)
    ref_groups = Dict{Tuple{String,String,Int}, Vector{Int}}()
    for i in eachindex(ref_v)
        push!(get!(ref_groups, (ref_v[i], ref_j[i], ref_len[i]), Int[]), i)
    end

    q_cdr3 = String.(coalesce.(query_table[!, cdr3_col], ""))
    q_v    = String.(coalesce.(query_table.v_call, ""))
    q_j    = String.(coalesce.(query_table.j_call, ""))
    q_ids  = hasproperty(query_table, :sequence_id) ?
        String.(coalesce.(query_table.sequence_id, "")) :
        ["query_$i" for i in 1:nrow(query_table)]

    results = ClonoQueryResult[]
    for qi in 1:nrow(query_table)
        isempty(q_v[qi]) && continue
        indices = get(ref_groups, (q_v[qi], q_j[qi], length(q_cdr3[qi])), Int[])
        matched = filter(ri -> is_similar_with_junction(
            q_cdr3[qi], ref_cdr3[ri], cq.mismatches, cq.cdr3_core), indices)
        match_df = isempty(matched) ? ref[1:0, :] : ref[matched, :]
        push!(results, ClonoQueryResult(q_ids[qi], length(matched), match_df))
    end
    results
end

# ─── Summary statistics ───

const CLONOQUERY_SUMMARY_COLS = [
    :FR1_SHM, :CDR1_SHM, :FR2_SHM, :CDR2_SHM, :FR3_SHM, :V_SHM, :J_SHM,
    :V_aa_mut, :J_aa_mut,
]

"""
    clonoquery_summary(results, ref_table) -> DataFrame

Average SHM values across matches for each query.
"""
function clonoquery_summary(results::Vector{ClonoQueryResult}, ref_table::DataFrame)
    available = filter(c -> hasproperty(ref_table, c), CLONOQUERY_SUMMARY_COLS)
    rows = NamedTuple[]
    for r in results
        avgs = Dict{Symbol,Float64}()
        for col in available
            if nrow(r.matches) > 0 && hasproperty(r.matches, col)
                vals = collect(skipmissing(r.matches[!, col]))
                avgs[col] = isempty(vals) ? 0.0 : sum(vals) / length(vals)
            else
                avgs[col] = 0.0
            end
        end
        push!(rows, (name=r.query_id, size=r.match_count,
                     NamedTuple{Tuple(Symbol.("avg_" .* string.(available))...)}(
                         Tuple(round(avgs[c]; digits=2) for c in available))...))
    end
    isempty(rows) ? DataFrame() : DataFrame(rows)
end

"""
    write_clonoquery(io, results)

Write results in reference format: comment lines with query info followed by matches.
"""
function write_clonoquery(io::IO, results::Vector{ClonoQueryResult})
    first_result = true
    for r in results
        first_result || println(io)
        first_result = false
        println(io, "# Query: $(r.query_id)")
        if nrow(r.matches) > 0
            for ri in 1:nrow(r.matches)
                println(io, join([string(r.matches[ri, c]) for c in names(r.matches)], '\t'))
            end
        end
    end
end
