# Configuration types and loading

struct PreprocessingFilter
    v_coverage::Float64
    j_coverage::Float64
    v_evalue::Float64
end

PreprocessingFilter(d::Dict) = PreprocessingFilter(
    Float64(d["v_coverage"]),
    Float64(d["j_coverage"]),
    Float64(d["v_evalue"]),
)

struct GermlineFilterCriteria
    unique_cdr3s::Int
    unique_js::Int
    whitelist::Bool
    cluster_size::Int
    allow_stop::Bool
    cross_mapping_ratio::Float64
    clonotype_ratio::Float64
    exact_ratio::Float64
    cdr3_shared_ratio::Float64
    unique_d_ratio::Float64
    unique_d_threshold::Int
end

GermlineFilterCriteria(d::Dict) = GermlineFilterCriteria(
    Int(d["unique_cdr3s"]),
    Int(d["unique_js"]),
    Bool(d["whitelist"]),
    Int(d["cluster_size"]),
    Bool(d["allow_stop"]),
    Float64(d["cross_mapping_ratio"]),
    Float64(d["clonotype_ratio"]),
    Float64(d["exact_ratio"]),
    Float64(d["cdr3_shared_ratio"]),
    Float64(d["unique_d_ratio"]),
    Int(d["unique_d_threshold"]),
)

struct JDiscoveryConfig
    allele_ratio::Float64
    cross_mapping_ratio::Float64
    propagate::Bool
end

JDiscoveryConfig(d::Dict) = JDiscoveryConfig(
    Float64(d["allele_ratio"]),
    Float64(d["cross_mapping_ratio"]),
    Bool(d["propagate"]),
)

"""
    Config

Complete IgDiscover configuration. Loaded from a TOML file merged with defaults.
"""
struct Config
    species::String
    sequence_type::String
    merge_program::String
    flash_maximum_overlap::Int
    limit::Int
    multialign_program::String
    minimum_merged_read_length::Int
    mismatch_penalty::Int
    barcode_length::Int
    barcode_consensus::Bool
    group_trim_g::Bool
    group_by_cdr3::String
    pseudo_cdr3_start::Int
    pseudo_cdr3_end::Int
    group_minimum_length::Int
    iterations::Int
    ignore_j::Bool
    d_coverage::Int
    subsample::Int
    consensus_threshold::Float64
    stranded::Bool
    rename::Bool
    race_g::Bool
    seed::Int
    exact_copies::Int
    preprocessing_filter::PreprocessingFilter
    pre_germline_filter::GermlineFilterCriteria
    germline_filter::GermlineFilterCriteria
    j_discovery::JDiscoveryConfig
end

"""
    load_config(path) -> Config

Load configuration from a TOML file, merging with built-in defaults.
"""
function load_config(path::AbstractString)
    defaults_path = joinpath(@__DIR__, "..", "config", "defaults.toml")
    defaults = TOML.parsefile(defaults_path)
    user = TOML.parsefile(path)
    parse_config(merge_config(defaults, user))
end

function load_config()
    path = "igdiscover.toml"
    isfile(path) || error("Configuration file '$path' not found. Run init_analysis first.")
    load_config(path)
end

function merge_config(defaults::Dict, user::Dict)
    result = deepcopy(defaults)
    for (k, v) in user
        if v isa Dict && haskey(result, k) && result[k] isa Dict
            merge!(result[k], v)
        else
            result[k] = v
        end
    end
    result
end

function parse_config(d::Dict)
    Config(
        get(d, "species", ""),
        get(d, "sequence_type", "Ig"),
        get(d, "merge_program", "pear"),
        Int(get(d, "flash_maximum_overlap", 300)),
        Int(get(d, "limit", 0)),
        get(d, "multialign_program", "muscle-fast"),
        Int(get(d, "minimum_merged_read_length", 300)),
        Int(get(d, "mismatch_penalty", 0)),
        Int(get(d, "barcode_length", 0)),
        Bool(get(d, "barcode_consensus", true)),
        Bool(get(d, "group_trim_g", false)),
        get(d, "group_by_cdr3", "none"),
        Int(get(d, "pseudo_cdr3_start", 80)),
        Int(get(d, "pseudo_cdr3_end", 60)),
        Int(get(d, "group_minimum_length", 0)),
        Int(get(d, "iterations", 1)),
        Bool(get(d, "ignore_j", false)),
        Int(get(d, "d_coverage", 70)),
        Int(get(d, "subsample", 1000)),
        Float64(get(d, "consensus_threshold", 60.0)),
        Bool(get(d, "stranded", false)),
        Bool(get(d, "rename", false)),
        Bool(get(d, "race_g", false)),
        Int(get(d, "seed", 1)),
        Int(get(d, "exact_copies", 0)),
        PreprocessingFilter(d["preprocessing_filter"]),
        GermlineFilterCriteria(d["pre_germline_filter"]),
        GermlineFilterCriteria(d["germline_filter"]),
        JDiscoveryConfig(d["j_discovery"]),
    )
end

"""
    write_default_config(path)

Write a default configuration file to the given path.
"""
function write_default_config(path::AbstractString)
    defaults_path = joinpath(@__DIR__, "..", "config", "defaults.toml")
    cp(defaults_path, path; force=false)
end
