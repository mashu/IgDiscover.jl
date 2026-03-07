module IgDiscover

using CSV
using CodecZlib
using Clustering
using DataFrames
using FASTX
using JSON3
using PrecompileTools
using Printf
using Random
using SHA
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
include("pipeline.jl")

export Config, load_config, write_default_config
export PreprocessingFilter, GermlineFilterCriteria, JDiscoveryConfig
export FastaRecord, read_fasta, read_fasta_dict, write_fasta, write_fasta_gz
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
export ClonotypeCaller, call_clonotypes

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
        # Shared utilities
        tallies(["a", "b", "a"])
        most_common(["a", "b", "a"])
        safe_divide(10, 5)
        safe_divide(1, 0)
        gene_family("IGHV1-2*01")
        is_same_gene("IGHV1-2*01", "IGHV1-2*02")
        from_toml(PreprocessingFilter, Dict("v_coverage" => 90.0, "j_coverage" => 60.0, "v_evalue" => 1e-3))

        # DataFrame column helpers
        df_test = DataFrame(a=["x", missing], b=[1, missing])
        ensure_column!(df_test, :a, "")
        ensure_column!(df_test, :b, 0)
        ensure_column!(df_test, :c, 0.0)

        # DNA utilities
        edit_distance("ATCGATCG", "ATCGATCA")
        edit_distance("ATCGATCG", "ATCGATCA"; maxdiff=2)
        hamming_distance("ATCG", "ATCA")
        translate("ATGATGATG")
        reverse_complement("ATCG")
        has_stop("ATGATG")
        sequence_hash("ATCGATCG")

        # FASTA I/O
        recs = read_fasta(fasta_path)
        d = read_fasta_dict(fasta_path)
        out_path = joinpath(tmpdir, "out.fasta")
        write_fasta(out_path, recs)

        # Config
        if isfile(defaults_path)
            raw = TOML.parsefile(defaults_path)
            parse_config(raw)
        end

        # Consensus
        consensus_sequence(["ATCG", "ATCG", "ATCA"]; threshold=0.6)

        # Header parsing
        parse_header("read1;size=5;barcode=ACG;")

        # CDR3
        cdr3_start_in_v("AAA" ^ 30 * "TTTTATTGTGCT", "IGH")
        cdr3_end_in_j("TGGGCAGGG", "IGH")

        # Alignment
        align_affine("ATCG", "ATCA")

        # Sanitization
        sanitize_imgt_sequence("ATG...CCC")
        allele_name_from_header("X|IGHV1*01|Homo")
    end
    rm(tmpdir; recursive=true, force=true)
end

end # module
