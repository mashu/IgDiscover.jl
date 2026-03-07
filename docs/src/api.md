# API Reference

## Core Types

```@docs
Config
FilterStats
```

[`FastaRecord`](@ref) and [`ClonotypeCaller`](@ref) are documented in the **FASTA I/O** and **Clonotypes** module pages.

## I/O

```@docs
read_fasta
read_fasta_dict
write_fasta
write_fasta_gz
read_assignments
write_table
write_table_gz
sanitize_imgt_record
sanitize_imgt_sequence
write_sanitized_imgt
```

## Pipeline

```@docs
init_analysis
run_pipeline
augment_table
filter_table
discover_germline
germline_filter!
deduplicate_by_consensus
germline_filter_to_fasta
discover_j_genes
rename_genes
call_clonotypes
```

For DNA utilities, alignment, CDR3 detection, and clustering see the module pages: **DNA Utilities**, **Alignment**, **CDR3 Detection**, **Clustering**.
