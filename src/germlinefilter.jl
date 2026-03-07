# Germline filter: filter V gene candidates to produce the final germline database.
# Uses abstract type + multiple dispatch for composable, extensible filters.

abstract type CandidateFilter end

# Immutable snapshot of a candidate used during pairwise filtering.
struct FilterCandidate
    sequence::String
    name::String
    clonotypes::Int
    exact::Int
    Ds_exact::Int
    cluster_size::Int
    whitelisted::Bool
    is_database::Bool
    cluster_size_is_accurate::Bool
    cdr3_start::Int
    index::Int
end

"""
    should_discard(filter, ref, candidate, same_gene) -> String

Return a non-empty reason if `candidate` should be discarded relative to `ref`.
Dispatches on the concrete filter type.
"""
function should_discard end

# ─── Filter implementations ───

struct IdenticalSequenceFilter <: CandidateFilter end

function should_discard(::IdenticalSequenceFilter, ref::FilterCandidate,
                        cand::FilterCandidate, ::Bool)
    cand.whitelisted && return ""
    ref.cluster_size < cand.cluster_size && return ""
    ref.sequence == cand.sequence      && return "identical_to=$(ref.name)"
    startswith(ref.sequence, cand.sequence) && return "identical_to=$(ref.name),truncated"
    ""
end

struct CrossMappingFilter <: CandidateFilter
    ratio::Float64
end

function should_discard(f::CrossMappingFilter, ref::FilterCandidate,
                        cand::FilterCandidate, ::Bool)
    s = SubString(ref.sequence,  1, min(ref.cdr3_start,  length(ref.sequence)))
    t = SubString(cand.sequence, 1, min(cand.cdr3_start, length(cand.sequence)))

    dist = if length(s) != length(t)
        short_len = min(length(s), length(t))
        min(
            edit_distance(s[1:short_len], t[1:short_len]; maxdiff=1),
            edit_distance(s[max(1,end-short_len+1):end], t[max(1,end-short_len+1):end]; maxdiff=1),
        )
    else
        edit_distance(s, t; maxdiff=1)
    end
    dist > 1 && return ""
    (!ref.is_database || !cand.is_database) && return ""

    total = ref.cluster_size + cand.cluster_size
    total == 0 && return ""
    ratio = cand.cluster_size / total
    (cand.cluster_size_is_accurate && ratio < f.ratio) ?
        @sprintf("xmap_ratio=%.4f,other=%s", ratio, ref.name) : ""
end

struct ClonotypeRatioFilter <: CandidateFilter
    ratio::Float64
end

function should_discard(f::ClonotypeRatioFilter, ref::FilterCandidate,
                        cand::FilterCandidate, same_gene::Bool)
    (!same_gene || ref.clonotypes == 0) && return ""
    ratio = cand.clonotypes / ref.clonotypes
    ratio < f.ratio ? @sprintf("clonotype_ratio=%.4f,other=%s", ratio, ref.name) : ""
end

struct ExactRatioFilter <: CandidateFilter
    ratio::Float64
end

function should_discard(f::ExactRatioFilter, ref::FilterCandidate,
                        cand::FilterCandidate, same_gene::Bool)
    (!same_gene || ref.exact == 0) && return ""
    ratio = cand.exact / ref.exact
    ratio < f.ratio ? @sprintf("ex_occ_ratio=%.1f,other=%s", ratio, ref.name) : ""
end

struct UniqueDRatioFilter <: CandidateFilter
    ratio::Float64
    threshold::Int
end

function should_discard(f::UniqueDRatioFilter, ref::FilterCandidate,
                        cand::FilterCandidate, same_gene::Bool)
    (!same_gene || ref.cluster_size < cand.cluster_size || ref.Ds_exact < f.threshold) && return ""
    ref.Ds_exact == 0 && return ""
    ratio = cand.Ds_exact / ref.Ds_exact
    ratio < f.ratio ? @sprintf("Ds_exact_ratio=%.1f,other=%s", ratio, ref.name) : ""
end

# ─── Whitelist ───

struct Whitelist
    sequences::Dict{String,String}   # sequence → name
end

Whitelist() = Whitelist(Dict{String,String}())

function Whitelist(paths::Vector{String})
    seqs = Dict{String,String}()
    for path in paths
        isfile(path) || continue
        for r in read_fasta(path)
            seqs[r.sequence] = r.name
        end
    end
    Whitelist(seqs)
end

# ─── Main filter ───

"""
    germline_filter!(candidates, criteria; whitelist) -> (filtered_df, annotated_df)

Apply per-entry and pairwise germline filters to V gene candidates.
"""
function germline_filter!(
    candidates::DataFrame,
    criteria::GermlineFilterCriteria;
    whitelist::Whitelist=Whitelist(),
)
    n = nrow(candidates)
    candidates.is_filtered  = zeros(Int, n)
    candidates.why_filtered = fill("", n)
    candidates.whitelist_diff = [haskey(whitelist.sequences, s) ? 0 : -1
                                  for s in candidates.consensus]

    mark_filtered!(candidates, candidates.CDR3s_exact .< criteria.unique_cdr3s,
                   "too_low_CDR3s_exact")
    mark_filtered!(candidates, candidates.Js_exact .< criteria.unique_js,
                   "too_low_Js_exact")

    if hasproperty(candidates, :CDR3_shared_ratio)
        ratios = [parse(Float64, string(x)) for x in candidates.CDR3_shared_ratio]
        mark_filtered!(candidates, ratios .> criteria.cdr3_shared_ratio,
                       "too_high_CDR3_shared_ratio")
    end

    if !criteria.allow_stop
        mark_filtered!(candidates,
            (candidates.has_stop .!= 0) .& (candidates.whitelist_diff .!= 0), "has_stop")
    end

    mark_filtered!(candidates,
        (candidates.cluster_size .< criteria.cluster_size) .& (candidates.whitelist_diff .!= 0),
        "too_low_cluster_size")

    @info "$(sum(candidates.is_filtered .== 0)) remain after per-entry filtering"

    filters = CandidateFilter[IdenticalSequenceFilter()]
    criteria.cross_mapping_ratio > 0 && push!(filters, CrossMappingFilter(criteria.cross_mapping_ratio))
    criteria.clonotype_ratio     > 0 && push!(filters, ClonotypeRatioFilter(criteria.clonotype_ratio))
    criteria.exact_ratio         > 0 && push!(filters, ExactRatioFilter(criteria.exact_ratio))
    criteria.unique_d_ratio      > 0 && push!(filters, UniqueDRatioFilter(criteria.unique_d_ratio, criteria.unique_d_threshold))

    has_db_diff = hasproperty(candidates, :database_diff)
    has_cdr3_start = hasproperty(candidates, :CDR3_start)

    fc = [FilterCandidate(
        candidates.consensus[i], candidates.name[i],
        candidates.clonotypes[i], candidates.exact[i], candidates.Ds_exact[i],
        candidates.cluster_size[i],
        candidates.whitelist_diff[i] == 0,
        has_db_diff && !ismissing(candidates.database_diff[i]) &&
            candidates.database_diff[i] == 0,
        occursin("all", candidates.cluster[i]) || occursin("db", candidates.cluster[i]),
        has_cdr3_start ? candidates.CDR3_start[i] : 10000,
        i,
    ) for i in 1:n]

    for ref in fc
        candidates.is_filtered[ref.index] > 0 && continue
        for cand in fc
            ref.index == cand.index && continue
            sg = is_same_gene(ref.name, cand.name)
            for f in filters
                reason = should_discard(f, ref, cand, sg)
                if !isempty(reason)
                    candidates.why_filtered[cand.index] *= reason * ";"
                    candidates.is_filtered[cand.index]  += 1
                    break
                end
            end
        end
    end

    annotated = copy(candidates)
    filtered  = candidates[candidates.is_filtered .== 0, :]
    select!(filtered, Not([:is_filtered, :why_filtered]))

    @info "$(nrow(filtered)) sequences in new database"
    (filtered, annotated)
end

"""
    germline_filter_to_fasta(filtered, path)

Write filtered V gene candidates to a FASTA file.
"""
function germline_filter_to_fasta(filtered::DataFrame, path::AbstractString)
    write_fasta(path, [FastaRecord(row.name, row.consensus) for row in eachrow(filtered)])
end

# ─── Internal helpers ───

function mark_filtered!(df::DataFrame, mask::AbstractVector{Bool}, reason::String)
    for i in eachindex(mask)
        if mask[i]
            df.why_filtered[i] *= reason * ";"
            df.is_filtered[i]  += 1
        end
    end
    @info "Marked $(sum(mask)) candidates as '$reason'"
end
