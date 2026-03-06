# J gene candidate discovery
# Mirrors IgDiscover22 discoverjd step

const J_MIN_EXACT = 3
const J_ALLELE_SUFFIX_LEN = 4   # last N bases used to distinguish alleles

struct JCandidate
    name::String
    source::String
    exact::Int
    CDR3s_exact::Int
    allele_ratio::Float64
    cross_mapping_ratio::Float64
    consensus::String
end

"""
    discover_j_genes(table, j_database, config) -> DataFrame

Discover J gene candidates from a filtered assignment table.
Groups reads by j_call, computes consensus for those with J_errors == 0,
applies allele-ratio and cross-mapping filters.
Returns a DataFrame with columns: name, source, exact, CDR3s_exact, consensus.
"""
function discover_j_genes(
    table::DataFrame,
    j_database::Dict{String,String},
    config::Config,
)
    hasproperty(table, :j_call)  || error("Table missing j_call column")
    hasproperty(table, :J_errors) || error("Table missing J_errors column")

    # Ensure required columns with defaults
    for col in (:cdr3,)
        hasproperty(table, col) || (table[!, col] = fill("", nrow(table)))
        table[!, col] = String.(coalesce.(table[!, col], ""))
    end

    candidates = JCandidate[]
    namer = NameGenerator()

    for gdf in groupby(table, :j_call)
        gene = String(first(gdf.j_call))
        isempty(gene) && continue

        exact_mask  = gdf.J_errors .== 0
        exact_group = gdf[exact_mask, :]
        nrow(exact_group) < J_MIN_EXACT && continue

        # Consensus from J_nt (ungapped J sequence alignment)
        j_seqs = if hasproperty(gdf, :j_sequence_alignment)
            filter(!isempty, String.(replace.(
                coalesce.(exact_group.j_sequence_alignment, ""), "-" => "")))
        else
            String[]
        end
        isempty(j_seqs) && continue

        cons = iterative_consensus(j_seqs;
            program="muscle-medium",
            threshold=config.consensus_threshold / 100,
            maximum_subsample_size=config.subsample)
        isempty(cons) && continue

        db_seq = get(j_database, gene, "")
        # If consensus matches database exactly (or is a prefix/suffix), keep db name
        seq_id = if !isempty(db_seq) && (cons == db_seq ||
                    startswith(db_seq, cons) || startswith(cons, db_seq))
            gene
        else
            unique_name(gene, cons)
        end

        unique_cdr3s = length(unique(filter(!isempty, String.(exact_group.cdr3))))

        push!(candidates, JCandidate(
            namer(seq_id), gene,
            nrow(exact_group), unique_cdr3s,
            0.0, 0.0, cons))
    end

    isempty(candidates) && return DataFrame()

    # Allele ratio filter: for alleles of the same gene, discard those with
    # exact count ratio below threshold relative to the best allele
    candidates = filter_j_alleles(candidates, config.j_discovery.allele_ratio)

    # Cross-mapping filter
    candidates = filter_j_cross_mapping(candidates, config.j_discovery.cross_mapping_ratio)

    @info "J gene discovery: $(length(candidates)) candidates"
    j_candidates_to_dataframe(candidates)
end

"""
    filter_j_alleles(candidates, allele_ratio) -> Vector{JCandidate}

Remove alleles whose exact-count is less than `allele_ratio` times the best allele's count
for the same gene family.
"""
function filter_j_alleles(candidates::Vector{JCandidate}, allele_ratio::Float64)
    # Group by gene family (part before '*')
    best_exact = Dict{String,Int}()
    for c in candidates
        family = first(split(c.source, '*'))
        best_exact[family] = max(get(best_exact, family, 0), c.exact)
    end
    filter(candidates) do c
        family = first(split(c.source, '*'))
        best = get(best_exact, family, 0)
        best == 0 || c.exact / best >= allele_ratio
    end
end

"""
    filter_j_cross_mapping(candidates, cross_mapping_ratio) -> Vector{JCandidate}

Remove candidates whose consensus is nearly identical to a more-expressed candidate
(edit distance ≤ 1 in the last J_ALLELE_SUFFIX_LEN bases).
"""
function filter_j_cross_mapping(candidates::Vector{JCandidate}, ratio::Float64)
    sorted = sort(candidates; by=c -> c.exact, rev=true)
    keep = trues(length(sorted))
    for i in eachindex(sorted)
        keep[i] || continue
        for j in (i+1):length(sorted)
            keep[j] || continue
            total = sorted[i].exact + sorted[j].exact
            total == 0 && continue
            # Compare suffix of consensus sequences
            si = sorted[i].consensus
            sj = sorted[j].consensus
            si_suf = length(si) >= J_ALLELE_SUFFIX_LEN ? si[end-J_ALLELE_SUFFIX_LEN+1:end] : si
            sj_suf = length(sj) >= J_ALLELE_SUFFIX_LEN ? sj[end-J_ALLELE_SUFFIX_LEN+1:end] : sj
            length(si_suf) == length(sj_suf) || continue
            hamming_distance(si_suf, sj_suf) <= 1 || continue
            xmap_ratio = sorted[j].exact / total
            xmap_ratio < ratio && (keep[j] = false)
        end
    end
    sorted[keep]
end

function j_candidates_to_dataframe(cs::Vector{JCandidate})
    DataFrame(
        name       = [c.name for c in cs],
        source     = [c.source for c in cs],
        exact      = [c.exact for c in cs],
        CDR3s_exact = [c.CDR3s_exact for c in cs],
        consensus  = [c.consensus for c in cs],
    )
end

"""
    discover_j_to_fasta(table, j_database, config, output_path) -> String

Run J discovery and write the resulting sequences to a FASTA file.
Returns output_path.
"""
function discover_j_to_fasta(
    table::DataFrame,
    j_database::Dict{String,String},
    config::Config,
    output_path::AbstractString,
)
    df = discover_j_genes(table, j_database, config)
    records = if nrow(df) == 0
        FastaRecord[]
    else
        [FastaRecord(row.name, row.consensus) for row in eachrow(df)]
    end
    write_fasta(output_path, records)
    output_path
end
