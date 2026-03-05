module IgDiscover

using CSV
using CodecZlib
using Clustering
using DataFrames
using Dates
using FASTX
using JSON3
using Logging
using Printf
using Random
using SHA
using Statistics
using TOML
using TranscodingStreams

include("dna.jl")
include("config.jl")
include("io.jl")
include("cdr3.jl")
include("alignment.jl")
include("clustering.jl")
include("igblast.jl")
include("augment.jl")
include("tablefilter.jl")
include("discovery.jl")
include("germlinefilter.jl")
include("pipeline.jl")

export Config, load_config, write_default_config
export PreprocessingFilter, GermlineFilterCriteria
export FastaRecord, read_fasta, read_fasta_dict, write_fasta
export read_assignments, write_table, write_table_gzipped
export multialign, consensus_sequence, iterative_consensus
export run_igblast_on_fasta, make_blastdb
export augment_table
export filter_table
export discover_germline
export germline_filter!, germline_filter_to_fasta
export init_analysis, run_pipeline

end # module
