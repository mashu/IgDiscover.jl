# Clustering algorithms for V gene discovery

"""
    distance_matrix(sequences; band=0.2) -> Matrix{Float64}

Pairwise edit distance matrix with bandwidth limit for efficiency.
Uses threading for the O(n²) unique-pair computation.
"""
function distance_matrix(sequences::Vector{String}; band::Float64 = 0.2)
    n = length(sequences)
    maxdiff = isempty(sequences) ? 0 : maximum(s -> round(Int, length(s) * band), sequences)

    # Deduplicate to avoid redundant comparisons
    unique_seqs = unique(sequences)
    nu = length(unique_seqs)

    # Compute pairwise distances for unique sequences (threaded)
    unique_dist_matrix = zeros(Int, nu, nu)
    pairs = [(i, j) for i in 1:nu for j in (i+1):nu]
    Threads.@threads for idx in eachindex(pairs)
        i, j = pairs[idx]
        d = min(maxdiff + 1, edit_distance(unique_seqs[i], unique_seqs[j]; maxdiff = maxdiff))
        unique_dist_matrix[i, j] = d
        unique_dist_matrix[j, i] = d
    end

    # Map unique sequences to indices for O(1) lookup
    seq_to_uid = Dict{String,Int}()
    for (i, s) in enumerate(unique_seqs)
        seq_to_uid[s] = i
    end

    M = zeros(Float64, n, n)
    @inbounds for i in 1:n, j in (i+1):n
        ui = seq_to_uid[sequences[i]]
        uj = seq_to_uid[sequences[j]]
        d = ui == uj ? 0.0 : Float64(unique_dist_matrix[ui, uj])
        M[i, j] = M[j, i] = d
    end
    M
end

"""
    cluster_sequences(sequences; minsize=5) -> (Matrix{Float64}, Vector{Int})

Hierarchical clustering with the igdiscover distance-ratio heuristic.
Returns (distance_matrix, cluster_ids) where cluster_id=0 means unassigned.
"""
function cluster_sequences(sequences::Vector{String}; minsize::Int = 5)
    n = length(sequences)
    n == 0 && return (zeros(0, 0), Int[])

    M = distance_matrix(sequences)
    hcl = Clustering.hclust(M; linkage = :average, branchorder = :optimal)
    clusters = zeros(Int, n)
    assign_igdiscover_clusters!(clusters, hcl, n, minsize)

    (M, clusters)
end

"""
Assign cluster IDs using igdiscover's distance-ratio heuristic.

Walk merge tree from highest to lowest distance. When ratio prev/current drops
below 0.8 and both subtrees have >= minsize leaves, assign a cluster.
"""
function assign_igdiscover_clusters!(clusters::Vector{Int},
                                    hcl::Clustering.Hclust,
                                    n::Int, minsize::Int)
    nmerges = length(hcl.heights)
    nmerges == 0 && return

    # Precompute leaf sets for each merge step
    node_leaves = Dict{Int,Vector{Int}}()

    function get_leaves(idx::Int)
        idx < 0 ? [-idx] : get(node_leaves, idx, Int[])
    end

    for step in 1:nmerges
        left, right = hcl.merges[step, 1], hcl.merges[step, 2]
        node_leaves[step] = vcat(get_leaves(left), get_leaves(right))
    end

    sorted_steps = sortperm(hcl.heights; rev = true)
    prev_dist = hcl.heights[sorted_steps[1]]
    cl = 1

    for step_idx in sorted_steps
        dist = hcl.heights[step_idx]
        left_leaves = get_leaves(hcl.merges[step_idx, 1])
        right_leaves = get_leaves(hcl.merges[step_idx, 2])

        if dist > 0 && prev_dist / dist < 0.8 &&
           length(left_leaves) >= minsize && length(right_leaves) >= minsize
            for id in left_leaves
                clusters[id] == 0 && (clusters[id] = cl)
            end
            cl += 1
        end
        prev_dist = dist
    end
end

# ─── Union-Find for single-linkage clustering ───

mutable struct UnionFind{T}
    parent::Dict{T,T}
    rank::Dict{T,Int}
end

UnionFind{T}() where T = UnionFind(Dict{T,T}(), Dict{T,Int}())

function add_node!(uf::UnionFind{T}, x::T) where T
    haskey(uf.parent, x) && return
    uf.parent[x] = x
    uf.rank[x] = 0
end

function find_root!(uf::UnionFind{T}, x::T) where T
    while uf.parent[x] != x
        uf.parent[x] = uf.parent[uf.parent[x]]  # path compression
        x = uf.parent[x]
    end
    x
end

function union!(uf::UnionFind{T}, x::T, y::T) where T
    rx, ry = find_root!(uf, x), find_root!(uf, y)
    rx == ry && return
    if uf.rank[rx] < uf.rank[ry]
        uf.parent[rx] = ry
    elseif uf.rank[rx] > uf.rank[ry]
        uf.parent[ry] = rx
    else
        uf.parent[ry] = rx
        uf.rank[rx] += 1
    end
end

function connected_components(uf::UnionFind{T}) where T
    groups = Dict{T,Vector{T}}()
    for x in keys(uf.parent)
        root = find_root!(uf, x)
        push!(get!(groups, root, T[]), x)
    end
    collect(values(groups))
end

"""
    single_linkage(strings, linked) -> Vector{Vector{String}}

Single-linkage clustering: `linked(s,t)` returns true if s and t should cluster together.
"""
function single_linkage(strings::Vector{String}, linked)
    uf = UnionFind{String}()
    for s in strings
        add_node!(uf, s)
    end
    for i in eachindex(strings), j in (i+1):length(strings)
        linked(strings[i], strings[j]) && union!(uf, strings[i], strings[j])
    end
    connected_components(uf)
end

"""
    count_clonotypes(cdr3s, j_calls; max_distance=6) -> Int

Cluster CDR3 sequences by J gene and edit distance, return number of clonotypes.
"""
function count_clonotypes(cdr3s::AbstractVector, j_calls::AbstractVector;
                         max_distance::Int = 6)
    total = 0
    j_groups = Dict{String,Vector{String}}()
    for (cdr3, j) in zip(cdr3s, j_calls)
        c = string(cdr3)
        isempty(c) && continue
        push!(get!(j_groups, string(j), String[]), c)
    end
    for group_cdr3s in values(j_groups)
        unique_cdr3s = unique(group_cdr3s)
        comps = single_linkage(unique_cdr3s,
            (s, t) -> edit_distance(s, t; maxdiff = max_distance) <= max_distance)
        total += length(comps)
    end
    total
end
