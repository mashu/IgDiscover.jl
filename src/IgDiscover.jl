module IgDiscover

using CSV
using CodecZlib
using Clustering
using DataFrames
using FASTX
using JSON3
using Printf
using Random
using SHA
using TOML
using TranscodingStreams

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

end # module
