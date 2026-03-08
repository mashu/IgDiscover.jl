module IgDiscover

import BioAlignments
import BioSequences
using CSV
using CodecZlib
using Clustering
using DataFrames
using FASTX
using JSON3
using PrecompileTools
using Printf
using ProgressMeter
using Random
using SHA
using ArgParse
using TOML
using TranscodingStreams

include("utils.jl")
include("dna.jl")
include("config.jl")
include("io.jl")
include("cdr3.jl")
include("alignment.jl")
include("clustering.jl")
include("group.jl")
include("igblast.jl")
include("augment.jl")
include("tablefilter.jl")
include("discovery.jl")
include("jdiscovery.jl")
include("germlinefilter.jl")
include("rename.jl")
include("clonotypes.jl")
include("clonoquery.jl")
include("plotalleles.jl")
include("upstream.jl")
include("count.jl")
include("dbdiff.jl")
include("haplotype.jl")
include("clusterplot.jl")
include("pipeline.jl")
include("cli.jl")

# ── Core pipeline ──
export Config, load_config, write_default_config
export PreprocessingFilter, GermlineFilterCriteria, JDiscoveryConfig
export FastaRecord, read_fasta, read_fasta_dict, read_fasta_dict_from_io
export has_fasta_records, write_fasta, write_fasta_gz
export sanitize_imgt_record, sanitize_imgt_sequence, write_sanitized_imgt
export read_assignments, write_table, write_table_gz
export multialign, consensus_sequence, iterative_consensus
export run_igblast_on_fasta, make_blastdb
export augment_table, parse_header
export filter_table, FilterStats
export discover_germline
export discover_j_genes, discover_j_to_fasta
export germline_filter!, germline_filter_to_fasta, deduplicate_by_consensus, Whitelist
export rename_genes
export init_analysis, run_pipeline
export group_reads

# ── Clonotypes ──
export ClonotypeCaller, call_clonotypes, write_members

# ── Clonoquery ──
export ClonoQuery, CDR3CoreSlice, ClonoQueryResult
export is_similar_with_junction, clonoquery_summary, write_clonoquery

# ── Allele usage ──
export AlleleUsage, AlleleFilterCriteria, allele_expression_counts

# ── Upstream analysis ──
export UpstreamParams, UpstreamAnalyzer

# ── Expression counts ──
export ExpressionCounter, natural_sort_key

# ── Database comparison ──
export DatabaseComparator, DatabaseDiff
export check_duplicate_names, check_duplicate_sequences, format_diff, exit_code

# ── Haplotype analysis ──
export HaplotypeAnalyzer, HaplotypePair, HaplotypeEntry
export format_tsv, switch_haplotype!, sort_haplotype!

# ── Cluster visualization ──
export ClusterPlotter, ClusterResult, render_heatmap

# ── Shared dispatch ──
export apply_filters, filter_by_allele_ratio

# ── Precompilation workload ──

@setup_workload begin
    tmpdir = mktempdir()
    fasta_path = joinpath(tmpdir, "test.fasta")
    open(fasta_path, "w") do io
        println(io, ">seq1\nATCGATCGATCG")
        println(io, ">seq2\nATCGATCGATCA")
        println(io, ">seq3\nGGCCGGCCGGCC")
    end
    defaults_path = joinpath(@__DIR__, "..", "config", "defaults.toml")

    @compile_workload begin
        tallies(["a", "b", "a"])
        most_common(["a", "b", "a"])
        safe_divide(10, 5)
        safe_divide(1, 0)
        gene_family("IGHV1-2*01")
        is_same_gene("IGHV1-2*01", "IGHV1-2*02")
        from_toml(PreprocessingFilter, Dict("v_coverage" => 90.0, "j_coverage" => 60.0, "v_evalue" => 1e-3))

        df_test = DataFrame(a=["x", missing], b=[1, missing])
        ensure_column!(df_test, :a, "")
        ensure_column!(df_test, :b, 0)
        ensure_column!(df_test, :c, 0.0)

        edit_distance("ATCGATCG", "ATCGATCA")
        edit_distance("ATCGATCG", "ATCGATCA"; maxdiff=2)
        hamming_distance("ATCG", "ATCA")
        translate("ATGATGATG")
        reverse_complement("ATCG")
        has_stop("ATGATG")
        sequence_hash("ATCGATCG")

        recs = read_fasta(fasta_path)
        d = read_fasta_dict(fasta_path)
        out_path = joinpath(tmpdir, "out.fasta")
        write_fasta(out_path, recs)

        if isfile(defaults_path)
            raw = TOML.parsefile(defaults_path)
            parse_config(raw)
        end

        consensus_sequence(["ATCG", "ATCG", "ATCA"]; threshold=0.6)
        parse_header("read1;size=5;barcode=ACG;")
        cdr3_start_in_v("AAA" ^ 30 * "TTTTATTGTGCT", "IGH")
        cdr3_end_in_j("TGGGCAGGG", "IGH")
        align_affine("ATCG", "ATCA")
        sanitize_imgt_sequence("ATG...CCC")
        allele_name_from_header("X|IGHV1*01|Homo")

        is_similar_with_junction("ATCGATCG", "ATCGATCA", 1.0, nothing)
        is_similar_with_junction("ATCGATCG", "ATCGATCA", 1.0, CDR3CoreSlice(2, 6))
        natural_sort_key("IGHV1-10*01")

        cmp = DatabaseComparator()
        cmp(recs[1:2], recs[1:2])
    end
    rm(tmpdir; recursive=true, force=true)
end

end # module
