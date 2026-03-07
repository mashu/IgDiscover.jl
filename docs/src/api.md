# API Reference

## Core Types

```@docs
FastaRecord
Config
FilterStats
Candidate
ClonotypeCaller
```

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

## DNA Utilities

```@docs
edit_distance
hamming_distance
translate
reverse_complement
has_stop
sequence_hash
```

## Alignment

```@docs
multialign
consensus_sequence
iterative_consensus
align_affine
describe_nt_change
```

## CDR3

```@docs
cdr3_start_in_v
cdr3_end_in_j
find_cdr3
```

## Clustering

```@docs
cluster_sequences
single_linkage
count_clonotypes
```
