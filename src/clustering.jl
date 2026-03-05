# Clustering algorithms for V gene discovery

"""
    distance_matrix(sequences; band=0.2) -> Matrix{Float64}

Pairwise edit distance matrix with bandwidth limit for efficiency.
"""
function distance_matrix(sequences::Vector{String}; band::Float64 = 0.2)
    n = length(sequences)
    maxdiff = isempty(sequences) ? 0 : maximum(s -> round(Int, length(s) * band), sequences)

    # Pre-compute distances between unique sequences
    unique_seqs = unique(sequences)
    unique_dists = Dict{Tuple{String,String},Int}()
    for i in eachindex(unique_seqs), j in (i+1):length(unique_seqs)
        s, t = unique_seqs[i], unique_seqs[j]
        d = min(maxdiff + 1, edit_distance(s, t; maxdiff = maxdiff))
        unique_dists[(s, t)] = d
        unique_dists[(t, s)] = d
    end

    M = zeros(Float64, n, n)
    for i in 1:n, j in (i+1):n
        s, t = sequences[i], sequences[j]
        d = s == t ? 0.0 : Float64(get(unique_dists, (s, t), maxdiff + 1))
        M[i, j] = M[j, i] = d
    end
    M
end

"""
    cluster_sequences(sequences; minsize=5) -> (Matrix{Float64}, Vector{Int})

Hierarchical clustering with the igdiscover distance-ratio heuristic.
Returns (distance_matrix, cluster_ids) where cluster_id=0 means unassigned.

The heuristic: walk inner nodes from highest to lowest distance.
When the ratio prev_dist/current_dist drops below 0.8, and both
subtrees have >= minsize leaves, assign the left subtree a cluster id.
"""
function cluster_sequences(sequences::Vector{String}; minsize::Int = 5)
    n = length(sequences)
    n == 0 && return (zeros(0, 0), Int[])

    M = distance_matrix(sequences)

    # Clustering.jl hclust: takes a full distance matrix
    hcl = Clustering.hclust(M; linkage = :average, branchorder = :optimal)

    clusters = zeros(Int, n)
    assign_igdiscover_clusters!(clusters, hcl, n, minsize)

    (M, clusters)
end

"""
Assign cluster IDs from hierarchical clustering using igdiscover's heuristic.

Clustering.jl's Hclust stores:
- merges: n-1 × 2 matrix. Negative values = leaf indices, positive = internal node indices.
- heights: n-1 vector of merge distances.
- order: leaf ordering.

We walk the merge tree from highest distance to lowest (like igdiscover),
splitting when the distance ratio drops below 0.8.
"""
function assign_igdiscover_clusters!(clusters::Vector{Int},
                                    hcl::Clustering.Hclust,
                                    n::Int, minsize::Int)
    nmerges = length(hcl.heights)
    nmerges == 0 && return

    # Build a mapping: for each merge step, which leaf indices belong to left and right
    node_leaves = Dict{Int,Vector{Int}}()

    function get_leaves(idx::Int)
        if idx < 0
            return [-idx]
        else
            return get(node_leaves, idx, Int[])
        end
    end

    for step in 1:nmerges
        left, right = hcl.merges[step, 1], hcl.merges[step, 2]
        node_leaves[step] = vcat(get_leaves(left), get_leaves(right))
    end

    # Sort merge steps by decreasing distance (highest first)
    sorted_steps = sortperm(hcl.heights; rev = true)

    prev_dist = hcl.heights[sorted_steps[1]]
    cl = 1

    for step_idx in sorted_steps
        dist = hcl.heights[step_idx]
        left, right = hcl.merges[step_idx, 1], hcl.merges[step_idx, 2]
        left_leaves = get_leaves(left)
        right_leaves = get_leaves(right)

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
