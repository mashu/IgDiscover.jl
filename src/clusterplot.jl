# Cluster plot — cluster sequences assigned to a gene and visualize
#
# Computes pairwise edit distances, performs hierarchical clustering,
# and produces a terminal heatmap via UnicodePlots or exports data for
# external visualization.

"""
    ClusterPlotter(; min_group_size=200, max_display=300, min_cluster_size=5,
                     gene_type=:V, ignore_j=false)

Callable struct that clusters sequences per gene and optionally renders
a Unicode heatmap to the terminal.

```julia
plotter = ClusterPlotter(min_group_size=100)
results = plotter(table)
for r in results
    println(r.gene, ": ", r.n_clusters, " clusters from ", r.n_sequences, " sequences")
end
```
"""
struct ClusterPlotter
    min_group_size::Int
    max_display::Int
    min_cluster_size::Int
    gene_type::Symbol
    ignore_j::Bool
end

ClusterPlotter(; min_group_size::Int=200, max_display::Int=300,
                 min_cluster_size::Int=5, gene_type::Symbol=:V,
                 ignore_j::Bool=false) =
    ClusterPlotter(min_group_size, max_display, min_cluster_size, gene_type, ignore_j)

struct ClusterResult
    gene::String
    n_sequences::Int
    n_clusters::Int
    cluster_ids::Vector{Int}
    distance_matrix::Matrix{Float64}
end

"""
    (plotter::ClusterPlotter)(table; genes=nothing) -> Vector{ClusterResult}

Cluster sequences for each gene. Optionally restrict to specific genes.
"""
function (plotter::ClusterPlotter)(table::DataFrame;
            genes::Union{Nothing,Vector{String}}=nothing)
    gene_col = gene_call_column(plotter.gene_type)
    seq_col = Symbol(string(plotter.gene_type), "_nt")

    df = copy(table)

    # Filter J errors unless ignored
    if !plotter.ignore_j && hasproperty(df, :J_SHM)
        df = df[coalesce.(df.J_SHM, Inf) .== 0.0, :]
    end

    hasproperty(df, gene_col) || error("Table missing column $gene_col")
    hasproperty(df, seq_col) || error("Table missing column $seq_col")

    results = ClusterResult[]
    gene_vals = String.(coalesce.(df[!, gene_col], ""))
    seq_vals = String.(coalesce.(df[!, seq_col], ""))

    gene_groups = Dict{String, Vector{Int}}()
    for i in eachindex(gene_vals)
        isempty(gene_vals[i]) && continue
        push!(get!(gene_groups, gene_vals[i], Int[]), i)
    end

    for gene in sort(collect(keys(gene_groups)))
        genes !== nothing && !(gene in genes) && continue
        indices = gene_groups[gene]
        length(indices) < plotter.min_group_size && continue

        sequences = seq_vals[indices]
        sequences = subsample(sequences, plotter.max_display)

        matrix = pairwise_distance_matrix(sequences)
        cluster_ids = hierarchical_cluster(matrix, plotter.min_cluster_size)

        push!(results, ClusterResult(gene, length(sequences),
              length(unique(filter(>(0), cluster_ids))), cluster_ids, matrix))
    end
    results
end

# ─── Subsampling (reservoir sampling) ───

function subsample(population::Vector{String}, size::Int)
    length(population) <= size && return copy(population)
    sample = population[1:size]
    for i in (size+1):length(population)
        j = rand(1:i)
        j <= size && (sample[j] = population[i])
    end
    sample
end

# ─── Pairwise distance matrix ───

function pairwise_distance_matrix(sequences::Vector{String}; band::Float64=0.2)
    n = length(sequences)
    maxdiff = maximum(s -> round(Int, length(s) * band), sequences; init=0)

    # Cache distances between unique sequences
    unique_seqs = unique(sequences)
    unique_dists = Dict{Tuple{String,String}, Float64}()
    for i in eachindex(unique_seqs), j in (i+1):length(unique_seqs)
        d = min(maxdiff + 1, edit_distance(unique_seqs[i], unique_seqs[j]; maxdiff=maxdiff))
        unique_dists[(unique_seqs[i], unique_seqs[j])] = Float64(d)
        unique_dists[(unique_seqs[j], unique_seqs[i])] = Float64(d)
    end

    matrix = zeros(Float64, n, n)
    for i in 1:n, j in (i+1):n
        d = sequences[i] == sequences[j] ? 0.0 :
            get(unique_dists, (sequences[i], sequences[j]), Float64(maxdiff + 1))
        matrix[i, j] = d
        matrix[j, i] = d
    end
    matrix
end

# ─── Hierarchical clustering with gap detection ───

function hierarchical_cluster(matrix::Matrix{Float64}, min_cluster_size::Int)
    n = size(matrix, 1)
    n < 2 && return zeros(Int, n)

    # Convert to condensed form for Clustering.jl
    dists = Float64[]
    for i in 1:n, j in (i+1):n
        push!(dists, matrix[i, j])
    end

    result = Clustering.hclust(matrix; linkage=:average, branchorder=:optimal)
    clusters = zeros(Int, n)

    # Walk the merge history and detect large distance jumps
    heights = result.heights
    isempty(heights) && return clusters

    max_height = maximum(heights)
    cl = 1
    prev_height = max_height

    # Build tree membership for gap-based splitting
    for (k, h) in enumerate(sort(heights; rev=true))
        h > 0 || continue
        if prev_height / h < 0.8
            # Large gap detected — assign a cluster
            members = Clustering.cutree(result; k=k)
            for group_id in unique(members)
                group_members = findall(==(group_id), members)
                if length(group_members) >= min_cluster_size &&
                   all(==(0), clusters[group_members])
                    for idx in group_members
                        clusters[idx] = cl
                    end
                    cl += 1
                end
            end
        end
        prev_height = h
    end

    clusters
end

# ─── Terminal visualization ───

"""
    render_heatmap(result::ClusterResult; width=60, height=30) -> String

Render a Unicode heatmap of the distance matrix. Falls back to a text
summary if UnicodePlots is not available.
"""
function render_heatmap(result::ClusterResult; width::Int=60, height::Int=30)
    n = size(result.distance_matrix, 1)
    lines = String[]
    push!(lines, "$(result.gene): $(result.n_sequences) sequences, $(result.n_clusters) clusters")
    push!(lines, "Distance matrix $(n)×$(n)")

    # Simple ASCII heatmap fallback
    step = max(1, n ÷ height)
    chars = [' ', '░', '▒', '▓', '█']
    max_d = maximum(result.distance_matrix)
    max_d == 0 && (max_d = 1.0)

    for i in 1:step:n
        row = Char[]
        col_step = max(1, n ÷ width)
        for j in 1:col_step:n
            d = result.distance_matrix[i, j]
            level = clamp(round(Int, d / max_d * 4) + 1, 1, 5)
            push!(row, chars[level])
        end
        push!(lines, String(row))
    end
    join(lines, '\n')
end
