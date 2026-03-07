# V gene candidate discovery — core algorithm

const MINGROUPSIZE = 5
const MINEXPRESSED = 10
const MAX_SUBSAMPLE = 1600

struct Candidate
    name::String
    source::String
    chain::String
    cluster::String
    cluster_size::Int
    Js::Int
    CDR3s::Int
    exact::Int
    full_exact::Int
    barcodes_exact::Int
    Ds_exact::Int
    Js_exact::Int
    CDR3s_exact::Int
    clonotypes::Int
    CDR3_exact_ratio::Float64
    CDR3_shared_ratio::Float64
    N_bases::Int
    database_diff::Int       # -1 = not in database
    database_changes::String
    has_stop::Int
    CDR3_start::Int
    consensus::String
end

struct DiscoveryParams
    database::Dict{String,String}
    windows::Vector{Tuple{Float64,Float64}}
    do_cluster::Bool
    consensus_threshold::Float64
    downsample::Int
    clonotype_differences::Int
    cluster_subsample_size::Int
    max_n_bases::Int
    exact_copies::Int
    d_coverage::Int
    d_evalue::Float64
    seed::Int
    cdr3_counts::Dict{String,Int}
end

struct GroupStats
    count::Int
    unique_D::Int
    unique_J::Int
    unique_CDR3::Int
    shared_CDR3_ratio::Float64
    clonotypes::Int
    unique_barcodes::Int
end

# ─── Name generator (functor — accumulates state) ───

mutable struct NameGenerator
    seen::Set{String}
    NameGenerator() = new(Set{String}())
end

function (ng::NameGenerator)(name::String)
    candidate = name
    ext = 'A'
    while candidate in ng.seen
        ext > 'Z' && error("Too many duplicate names for '$name'")
        candidate = name * ext
        ext = Char(UInt32(ext) + 1)
    end
    push!(ng.seen, candidate)
    candidate
end

# ─── Helpers ───

function set_discovery_seed!(gene::String, base_seed::Int)
    h = bytes2hex(sha256(Vector{UInt8}(codeunits(gene))))[end-7:end]
    Random.seed!(parse(UInt32, h; base=16) + base_seed)
end

function compute_group_stats(
    group::DataFrame,
    params::DiscoveryParams,
    other_cdr3::Dict{String,Int},
)
    cdr3s = filter(!isempty, String.(group.cdr3))
    cdr3_counts = tallies(cdr3s)
    unique_cdr3 = length(cdr3_counts)
    shared = count(k -> get(other_cdr3, k, 0) > 0, keys(cdr3_counts))
    shared_ratio = unique_cdr3 == 0 ? 0.0 : shared / unique_cdr3

    unique_j = length(unique(filter(!isempty, String.(group.j_call))))

    clono = count_clonotypes(group.cdr3, group.j_call; max_distance=params.clonotype_differences)

    d_mask = (group.D_errors .== 0) .&
             (coalesce.(group.D_covered, 0.0) .>= params.d_coverage) .&
             (coalesce.(group.d_support, Inf) .<= params.d_evalue)
    unique_d = length(unique(filter(!isempty, String.(group[d_mask, :d_call]))))

    unique_bc = hasproperty(group, :barcode) ?
        length(unique(filter(!isempty, String.(group.barcode)))) : 0

    GroupStats(nrow(group), unique_d, unique_j, unique_cdr3,
              shared_ratio, clono, unique_bc)
end

function sibling_consensus(gene::String, group::DataFrame, params::DiscoveryParams)
    seqs = filter(!isempty, String.(group.V_nt))
    isempty(seqs) && return ""
    cons = iterative_consensus(seqs;
        program="muscle-medium",
        threshold=params.consensus_threshold / 100,
        maximum_subsample_size=params.downsample)
    if haskey(params.database, gene)
        dbseq = params.database[gene]
        (startswith(cons, dbseq) || startswith(dbseq, cons)) && return dbseq
    end
    cons
end

"""
    merge_n_tolerant(s, t) -> Union{String, Nothing}

Merge two sequences where 'N' is treated as a wildcard. Returns nothing on conflict.
"""
function merge_n_tolerant(s::String, t::String)
    len = max(length(s), length(t))
    buf = Char[]
    for i in 1:len
        c1 = i <= length(s) ? s[i] : '-'
        c2 = i <= length(t) ? t[i] : '-'
        result = if c1 == '-'; c2
        elseif c2 == '-'; c1
        elseif c1 == 'N'; c2
        elseif c2 == 'N'; c1
        elseif c1 == c2; c1
        else return nothing
        end
        push!(buf, result)
    end
    String(buf)
end

# ─── Discovery for a single gene ───

function discover_gene(gene::String, assignments::DataFrame, params::DiscoveryParams)
    set_discovery_seed!(gene, params.seed)

    cdr3_vals = filter(>(0), assignments.V_CDR3_start)
    cdr3_start_val = isempty(cdr3_vals) ? 0 : most_common(cdr3_vals)[1]

    v_no_cdr3 = [cdr3_start_val > 0 && cdr3_start_val <= length(s) ?
                  s[1:cdr3_start_val] : s for s in assignments.V_nt]

    v_no_cdr3_index = Dict{String,Vector{Int}}()
    for (i, seq) in enumerate(v_no_cdr3)
        push!(get!(v_no_cdr3_index, seq, Int[]), i)
    end

    v_nt_index = Dict{String,Vector{Int}}()
    for (i, seq) in enumerate(assignments.V_nt)
        push!(get!(v_nt_index, seq, Int[]), i)
    end

    siblings = Tuple{String,String,Vector{Int}}[]

    for (left, right) in params.windows
        indices = findall(i -> left <= assignments.V_SHM[i] < right, 1:nrow(assignments))
        length(indices) < MINGROUPSIZE && continue
        group = assignments[indices, :]
        seq = sibling_consensus(gene, group, params)
        label = (left, right) == (0.0, 100.0) ? "all" :
            "$(left == floor(left) ? Int(left) : left)-$(right == floor(right) ? Int(right) : right)"
        push!(siblings, (seq, label, indices))
    end

    if params.do_cluster && nrow(assignments) >= MINGROUPSIZE
        sample_idx = reservoir_sample(collect(1:nrow(assignments)), params.cluster_subsample_size)
        _, cluster_ids = cluster_sequences(v_no_cdr3[sample_idx]; minsize=MINGROUPSIZE)

        groups_by_cl = Dict{Int,Vector{Int}}()
        for (i, cid) in enumerate(cluster_ids)
            cid == 0 && continue
            push!(get!(groups_by_cl, cid, Int[]), sample_idx[i])
        end

        cl_num = 0
        for idx_list in values(groups_by_cl)
            length(idx_list) < MINGROUPSIZE && continue
            seq = sibling_consensus(gene, assignments[idx_list, :], params)
            cl_num += 1
            push!(siblings, (seq, "cl$cl_num", idx_list))
        end
    end

    db_seq = get(params.database, gene, "")
    if !isempty(db_seq) && !any(s == db_seq for (s, _, _) in siblings)
        exact_idx = findall(==(0), assignments.V_errors)
        if length(exact_idx) >= MINEXPRESSED
            @debug "Re-adding expressed database sequence '$gene'"
            push!(siblings, (db_seq, "db", exact_idx))
        end
    end

    candidates = Candidate[]
    for (seq, sib_name, sib_indices) in siblings
        isempty(seq) && continue
        n_bases = count(==('N'), seq)
        n_bases > params.max_n_bases && continue

        v_no_cdr3_seq = cdr3_start_val > 0 && cdr3_start_val <= length(seq) ?
            seq[1:cdr3_start_val] : seq

        exact_idx    = get(v_no_cdr3_index, v_no_cdr3_seq, Int[])
        sib_group    = assignments[sib_indices, :]
        exact_group  = assignments[exact_idx, :]

        full_exact_idx = get(v_nt_index, seq, Int[])

        sib_cdr3   = tallies(filter(!isempty, String.(sib_group.cdr3)))
        other_cdr3 = copy(params.cdr3_counts)
        for (k, v) in sib_cdr3
            haskey(other_cdr3, k) && (other_cdr3[k] = max(0, other_cdr3[k] - v))
        end

        info_window = compute_group_stats(sib_group, params, other_cdr3)
        info_exact  = compute_group_stats(exact_group, params, other_cdr3)

        db_diff    = haskey(params.database, gene) ? edit_distance(seq, params.database[gene]) : -1
        db_changes = db_diff > 0 ? describe_nt_change(params.database[gene], seq) : ""

        db_diff > 0 && info_exact.count < 2 && continue

        seq_id     = db_diff == 0 ? gene : unique_name(gene, seq)
        chain_val  = isempty(sib_group.locus) ? "" : most_common(String.(sib_group.locus))[1]
        cdr3_ratio = safe_divide(info_exact.count, info_exact.unique_CDR3)

        push!(candidates, Candidate(
            seq_id, gene, chain_val, sib_name,
            info_window.count,
            info_window.unique_J, info_window.unique_CDR3,
            info_exact.count, length(full_exact_idx),
            info_exact.unique_barcodes,
            info_exact.unique_D, info_exact.unique_J, info_exact.unique_CDR3,
            info_exact.clonotypes, cdr3_ratio, info_exact.shared_CDR3_ratio,
            n_bases, db_diff, db_changes, has_stop(seq) ? 1 : 0,
            cdr3_start_val, seq))
    end
    candidates
end

# ─── Prepare table columns for discovery ───

function prepare_discovery_table!(table::DataFrame)
    for col in (:cdr3, :barcode, :j_call, :d_call, :locus)
        ensure_column!(table, col, "")
    end
    for col in (:cdr3, :j_call, :V_nt)
        table[!, col] = String.(coalesce.(table[!, col], ""))
    end
    for col in (:V_errors, :D_errors, :J_errors, :V_CDR3_start)
        ensure_column!(table, col, 0)
    end
    for col in (:V_SHM, :D_covered)
        ensure_column!(table, col, 0.0)
    end
    hasproperty(table, :d_support) || (table.d_support = fill(Inf, nrow(table)))
end

# ─── Main entry point ───

"""
    discover_germline(table, database, config) -> DataFrame

Discover V gene candidates from a filtered assignment table.
"""
function discover_germline(
    table::DataFrame,
    database::Dict{String,String},
    config::Config,
)
    table = copy(table)

    if !config.ignore_j && hasproperty(table, :J_SHM)
        table = table[table.J_SHM .== 0.0, :]
        @info "$(nrow(table)) rows remain after J%SHM=0 filter"
    end

    if !hasproperty(table, :V_nt) && hasproperty(table, :v_sequence_alignment)
        table.V_nt = replace.(coalesce.(table.v_sequence_alignment, ""), "-" => "")
    end
    hasproperty(table, :V_nt) || error("Table must have V_nt or v_sequence_alignment column")

    prepare_discovery_table!(table)

    cdr3_counts = tallies(filter(!isempty, String.(table.cdr3)))
    @info "$(length(cdr3_counts)) unique CDR3s overall"

    windows = [(Float64(s), Float64(s) + 2.0) for s in 0.0:2.0:18.0]
    push!(windows, (0.0, 100.0))

    params = DiscoveryParams(
        database, windows, true, config.consensus_threshold, MAX_SUBSAMPLE,
        6, config.subsample, 0, max(config.exact_copies, 1),
        config.d_coverage, 1e-4, config.seed, cdr3_counts)

    all_candidates = Candidate[]
    namer = NameGenerator()

    for gdf in groupby(table, :v_call)
        gene = String(first(gdf.v_call))
        isempty(gene) && continue
        nrow(gdf) < MINGROUPSIZE && continue
        for c in discover_gene(gene, DataFrame(gdf), params)
            named = Candidate(
                namer(c.name), c.source, c.chain, c.cluster,
                c.cluster_size, c.Js, c.CDR3s, c.exact, c.full_exact, c.barcodes_exact,
                c.Ds_exact, c.Js_exact, c.CDR3s_exact, c.clonotypes,
                c.CDR3_exact_ratio, c.CDR3_shared_ratio, c.N_bases,
                c.database_diff, c.database_changes, c.has_stop, c.CDR3_start, c.consensus)
            push!(all_candidates, named)
        end
    end

    @info "$(length(all_candidates)) candidate sequences generated"
    candidates_to_dataframe(all_candidates)
end

function candidates_to_dataframe(cs::Vector{Candidate})
    isempty(cs) && return DataFrame()
    DataFrame(
        name              = [c.name for c in cs],
        source            = [c.source for c in cs],
        chain             = [c.chain for c in cs],
        cluster           = [c.cluster for c in cs],
        cluster_size      = [c.cluster_size for c in cs],
        Js                = [c.Js for c in cs],
        CDR3s             = [c.CDR3s for c in cs],
        exact             = [c.exact for c in cs],
        full_exact        = [c.full_exact for c in cs],
        barcodes_exact    = [c.barcodes_exact for c in cs],
        Ds_exact          = [c.Ds_exact for c in cs],
        Js_exact          = [c.Js_exact for c in cs],
        CDR3s_exact       = [c.CDR3s_exact for c in cs],
        clonotypes        = [c.clonotypes for c in cs],
        CDR3_exact_ratio  = [@sprintf("%.2f", c.CDR3_exact_ratio) for c in cs],
        CDR3_shared_ratio = [@sprintf("%.2f", c.CDR3_shared_ratio) for c in cs],
        N_bases           = [c.N_bases for c in cs],
        database_diff     = [c.database_diff == -1 ? missing : c.database_diff for c in cs],
        database_changes  = [c.database_changes for c in cs],
        has_stop          = [c.has_stop for c in cs],
        CDR3_start        = [c.CDR3_start for c in cs],
        consensus         = [c.consensus for c in cs],
    )
end
