# Clonotype calling — groups sequences by V gene, J gene, and CDR3 similarity
#
# Two sequences belong to the same clonotype if they share:
#   - V gene assignment
#   - J gene assignment
#   - CDR3 nucleotide similarity (Hamming distance ≤ max_mismatches)
#
# Uses single-linkage clustering within each (V, J) group.

struct ClonotypeCaller
    max_mismatches::Int
    v_shm_threshold::Float64
    sort_by_size::Bool
end

ClonotypeCaller(; max_mismatches::Int=1, v_shm_threshold::Float64=5.0, sort_by_size::Bool=false) =
    ClonotypeCaller(max_mismatches, v_shm_threshold, sort_by_size)

struct ClonotypeMember
    row_index::Int
    cdr3_nt::String
    v_shm::Float64
end

"""
    call_clonotypes(table; caller) -> (clonotypes::DataFrame, members::Vector{Vector{Int}})

Assign clonotypes to a filtered assignment table.

Returns a DataFrame with one row per clonotype (the representative) plus a
`clonotype_id` column, and a vector of member row-index groups.
"""
function call_clonotypes(table::DataFrame;
                        caller::ClonotypeCaller=ClonotypeCaller())
    call_clonotypes(table, caller)
end

function call_clonotypes(table::DataFrame, caller::ClonotypeCaller)
    required = (:v_call, :j_call, :cdr3)
    for col in required
        hasproperty(table, col) || error("Table missing required column: $col")
    end

    # Ensure string columns
    v_calls = String.(coalesce.(table.v_call, ""))
    j_calls = String.(coalesce.(table.j_call, ""))
    cdr3s   = String.(coalesce.(table.cdr3, ""))
    v_shms  = hasproperty(table, :V_SHM) ?
        Float64.(coalesce.(table.V_SHM, 0.0)) :
        zeros(Float64, nrow(table))

    # Group by (V gene, J gene)
    vj_groups = Dict{Tuple{String,String}, Vector{ClonotypeMember}}()
    for i in 1:nrow(table)
        vc, jc, c3 = v_calls[i], j_calls[i], cdr3s[i]
        (isempty(vc) || isempty(jc) || isempty(c3)) && continue
        key = (vc, jc)
        push!(get!(vj_groups, key, ClonotypeMember[]),
              ClonotypeMember(i, c3, v_shms[i]))
    end

    all_clusters = Vector{Vector{Int}}()
    representatives = Int[]

    for ((vc, jc), members) in sort!(collect(vj_groups); by=first)
        clusters = cluster_cdr3s(members, caller.max_mismatches)
        for cluster_members in clusters
            push!(all_clusters, [m.row_index for m in cluster_members])
            rep = pick_representative(cluster_members)
            push!(representatives, rep.row_index)
        end
    end

    # Build output table
    clonotype_df = table[representatives, :]
    clonotype_df.clonotype_id = collect(1:length(representatives))
    clonotype_df.clonotype_size = [length(c) for c in all_clusters]

    # Compute member mutation rates if V_SHM available
    if hasproperty(table, :V_SHM)
        add_mutation_rates!(clonotype_df, table, all_clusters, caller.v_shm_threshold)
    end

    if caller.sort_by_size
        perm = sortperm(clonotype_df.clonotype_size; rev=true)
        clonotype_df = clonotype_df[perm, :]
        all_clusters = all_clusters[perm]
    end

    (clonotype_df, all_clusters)
end

"""
Cluster CDR3 sequences by Hamming distance using single-linkage.
Returns groups of ClonotypeMember.
"""
function cluster_cdr3s(members::Vector{ClonotypeMember}, max_mismatches::Int)
    isempty(members) && return Vector{Vector{ClonotypeMember}}()
    length(members) == 1 && return [members]

    # Group by CDR3 string for deduplication
    cdr3_to_members = Dict{String, Vector{ClonotypeMember}}()
    for m in members
        push!(get!(cdr3_to_members, m.cdr3_nt, ClonotypeMember[]), m)
    end

    unique_cdr3s = collect(keys(cdr3_to_members))

    if length(unique_cdr3s) == 1
        return [members]
    end

    # Single-linkage on unique CDR3s
    linked(s, t) = length(s) == length(t) && hamming_distance(s, t) <= max_mismatches
    cdr3_clusters = single_linkage(unique_cdr3s, linked)

    # Expand back to full member lists
    [reduce(vcat, cdr3_to_members[c] for c in cluster)
     for cluster in cdr3_clusters]
end

"""
Pick the least mutated member as clonotype representative.
"""
function pick_representative(members::Vector{ClonotypeMember})
    best = members[1]
    for m in members
        m.v_shm < best.v_shm && (best = m)
    end
    best
end

"""
Compute CDR3 and VDJ mutation rate columns for the clonotype table.
"""
function add_mutation_rates!(clonotype_df::DataFrame,
                            table::DataFrame,
                            clusters::Vector{Vector{Int}},
                            v_shm_threshold::Float64)
    n = nrow(clonotype_df)
    cdr3_nt_mindiffrate = fill(NaN, n)

    has_cdr3 = hasproperty(table, :cdr3)
    has_vshm = hasproperty(table, :V_SHM)

    for (ci, cluster_indices) in enumerate(clusters)
        length(cluster_indices) <= 1 && continue
        !has_cdr3 && continue

        # Find reference (lowest V_SHM in cluster)
        ref_idx = cluster_indices[1]
        ref_shm = has_vshm ? Float64(coalesce(table.V_SHM[ref_idx], 0.0)) : 0.0
        for idx in cluster_indices
            shm = has_vshm ? Float64(coalesce(table.V_SHM[idx], 0.0)) : 0.0
            if shm < ref_shm
                ref_shm = shm
                ref_idx = idx
            end
        end

        ref_shm > v_shm_threshold && continue

        ref_cdr3 = String(coalesce(table.cdr3[ref_idx], ""))
        isempty(ref_cdr3) && continue

        # For the representative row in clonotype_df, compute its own diffrate
        rep_cdr3 = String(coalesce(table.cdr3[cluster_indices[1]], ""))
        if !isempty(rep_cdr3) && !isempty(ref_cdr3)
            d = edit_distance(rep_cdr3, ref_cdr3)
            cdr3_nt_mindiffrate[ci] = 100.0 * d / length(ref_cdr3)
        end
    end

    clonotype_df.CDR3_nt_mindiffrate = cdr3_nt_mindiffrate
end

"""
    write_clonotypes(path, clonotype_df)

Write clonotype table to TSV.
"""
function write_clonotypes(path::AbstractString, df::DataFrame)
    write_table(path, df)
end

"""
    write_members(path, table, clusters)

Write members file: each clonotype's constituent rows, separated by blank lines.
"""
function write_members(path::AbstractString, table::DataFrame, clusters::Vector{Vector{Int}})
    open(path, "w") do io
        header = join(names(table), '\t')
        println(io, header)
        for (i, cluster_indices) in enumerate(clusters)
            i > 1 && println(io)
            for idx in cluster_indices
                row_vals = [string(table[idx, col]) for col in names(table)]
                println(io, join(row_vals, '\t'))
            end
        end
    end
end
