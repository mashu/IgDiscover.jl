# Output formats

This page describes every tabular (TSV/TSV.GZ) and JSON output produced by IgDiscover.jl: file names, column names, meaning, and how values are computed where applicable.

---

## Overview of output files

| File | Description |
|------|-------------|
| `airr.tsv.gz` | Raw IgBLAST assignment (AIRR format) |
| `assigned.tsv.gz` | AIRR + IgDiscover-added columns |
| `filtered.tsv.gz` | Subset of assigned after preprocessing filter |
| `candidates.tsv` | V gene discovery candidates (one row per candidate) |
| `new_V_pregermline.tsv` | V candidates after pre-germline filter |
| `new_V_germline.tsv` | V candidates after full germline filter (final V database source) |
| `annotated_V_pregermline.tsv` / `annotated_V_germline.tsv` | Same as candidates + filter annotations |
| `clonotypes.tsv` | One row per clonotype (representative sequence + metadata) |
| `stats/filtered.json` | Preprocessing filter counts |

FASTA outputs (`new_V_*.fasta`, `new_J.fasta`) contain one record per accepted sequence; headers are gene/allele names and the sequence is the consensus nucleotide.

---

## 1. `airr.tsv.gz`

IgBLAST output with **-outfmt 19** (AIRR Rearrangement schema). One row per input read.

Columns are defined by the [AIRR Community Rearrangement specification](https://docs.airr-community.org/). IgDiscover uses in particular:

| Column | Meaning |
|--------|---------|
| `sequence_id` | Query identifier (from FASTA header) |
| `v_call`, `d_call`, `j_call` | Assigned V/D/J gene (and allele) |
| `locus` | Chain type (e.g. IGH, IGK, IGL) |
| `v_identity`, `j_identity` | Percent identity of V/J alignment (0–100) |
| `v_germline_alignment`, `v_sequence_alignment` | Aligned reference and query strings (with gaps) |
| `d_germline_alignment`, `d_sequence_alignment` | Same for D |
| `j_germline_alignment`, `j_sequence_alignment` | Same for J |
| `v_germline_start`, `v_germline_end` | 1-based reference coordinates |
| `v_sequence_start`, `v_sequence_end` | 1-based query coordinates (and similarly for D, J) |
| `sequence` | Query nucleotide sequence |
| `stop_codon` | Whether translation has stop (e.g. true/false) |
| `d_support` | E-value for D hit (when present) |

IgBLAST may add other AIRR columns; only those used by IgDiscover are listed above.

---

## 2. `assigned.tsv.gz`

Same as `airr.tsv.gz` plus **IgDiscover-added columns** from [`augment_table`](@ref). One row per read.

### Columns added by IgDiscover

| Column | Type | Meaning / formula |
|--------|------|-------------------|
| `sequence_id` | string | **Overwritten** from header: first field before `;` (e.g. from `name;size=5;barcode=ACG` → `name`) |
| `count` | int | Group size from header: `size=N` → N; 0 if absent |
| `barcode` | string | From header: `barcode=...`; empty if absent |
| `V_SHM` | float | V gene somatic hypermutation: `V_SHM = 100 - v_identity` (percent difference from germline) |
| `J_SHM` | float | J gene somatic hypermutation: `J_SHM = 100 - j_identity` |
| `V_errors` | int | Number of position-wise differences in V alignment (mismatches + gaps): `count(g ≠ s for (g,s) in zip(v_germline_alignment, v_sequence_alignment))` |
| `D_errors` | int | Same for D alignment |
| `J_errors` | int | Same for J alignment |
| `V_covered` | float | Fraction of reference V length aligned (ungapped): `100 × (number of non-gap chars in v_germline_alignment) / length(reference_V)` |
| `D_covered` | float | Same for D |
| `J_covered` | float | Same for J |
| `V_CDR3_start` | int | 1-based position in the **query** sequence where the V-side CDR3 boundary falls (from locus-specific rules + alignment); 0 if not determined |
| `cdr3` | string | Nucleotide CDR3: `sequence[cdr3_start : cdr3_end]` when both boundaries are found |
| `cdr3_aa` | string | Amino acid CDR3: `translate(cdr3)` |
| `cdr3_start`, `cdr3_end` | int | 1-based query coordinates of CDR3 |
| `V_nt` | string | V region nucleotides (alignment gaps removed): `replace(v_sequence_alignment, "-" => "")` |
| `stop_codon` | string | Normalized to `"T"` / `"F"` (true/false) |
| `v_call`, `d_call`, `j_call` | string | **Overwritten** to short allele form (e.g. `IGHV1-18*01`) and leading `%` stripped |

All other columns from AIRR are passed through unchanged.

---

## 3. `filtered.tsv.gz`

Same **columns** as `assigned.tsv.gz`. Contains only rows that pass the preprocessing filter (see Configuration and Pipeline):

- Non-empty `v_call` and `j_call`
- No internal stop codon (`stop_codon` = false)
- V E-value ≤ `preprocessing_filter.v_evalue` (when `v_support` present)
- V coverage ≥ `preprocessing_filter.v_coverage` (%)
- J coverage ≥ `preprocessing_filter.j_coverage` (%)

No new columns are added.

---

## 4. `stats/filtered.json`

JSON object with one integer field per stage of the preprocessing filter:

| Field | Meaning |
|-------|---------|
| `total` | Total rows in the table before filtering |
| `has_vj_assignment` | Rows with both V and J assigned |
| `has_no_stop` | Rows after removing stop codons |
| `good_v_evalue` | Rows passing V E-value threshold |
| `good_v_coverage` | Rows passing V coverage threshold |
| `good_j_coverage` | Rows passing J coverage threshold |
| `has_cdr3` | Count of rows with non-empty CDR3 (informational) |

---

## 5. `candidates.tsv`

V gene **discovery** output: one row per candidate sequence. Columns come from [`discover_germline`](@ref) (see Discovery module).

| Column | Type | Meaning |
|--------|------|---------|
| `name` | string | Candidate identifier (gene/allele name or unique name with sequence suffix) |
| `source` | string | Gene name (e.g. IGHV1-18*01) |
| `chain` | string | Locus (e.g. IGH) |
| `cluster` | string | Cluster label (e.g. cl1, db) |
| `cluster_size` | int | Number of reads in the cluster used for this candidate |
| `Js` | int | Unique J gene assignments in the cluster |
| `CDR3s` | int | Unique CDR3 sequences in the cluster |
| `exact` | int | Reads with exact V match (no errors) in the info window |
| `full_exact` | int | Reads with exact full V sequence match |
| `barcodes_exact` | int | Unique barcodes among exact matches |
| `Ds_exact` | int | Unique D assignments among exact matches (with D coverage/evalue criteria) |
| `Js_exact` | int | Unique J among exact matches |
| `CDR3s_exact` | int | Unique CDR3 among exact matches |
| `clonotypes` | int | Clonotype count in the exact group (single-linkage by CDR3, max distance from config) |
| `CDR3_exact_ratio` | string | `exact_count / unique_CDR3` in exact group (formatted "%.2f") |
| `CDR3_shared_ratio` | string | Fraction of exact-group CDR3s that appear elsewhere in the dataset (formatted "%.2f") |
| `N_bases` | int | Number of N characters in the consensus sequence |
| `database_diff` | int | Edit distance from current database sequence for this gene; -1 or missing if not in database |
| `database_changes` | string | HGVS-style nucleotide change description vs database (e.g. 123A>G) |
| `has_stop` | int | 1 if consensus contains internal stop codon, else 0 |
| `CDR3_start` | int | 1-based position of CDR3 start in consensus (0 if N/A) |
| `consensus` | string | Consensus nucleotide sequence for this candidate |

---

## 6. `new_V_pregermline.tsv` and `new_V_germline.tsv`

Same **column set** as `candidates.tsv`, but only rows that pass the germline filter (per-entry and pairwise). The **pre** version uses `pre_germline_filter` criteria; the **germline** version uses `germline_filter` criteria. No extra columns; [`deduplicate_by_consensus`](@ref) keeps one row per unique `consensus` (largest `cluster_size` wins).

---

## 7. `annotated_V_pregermline.tsv` and `annotated_V_germline.tsv`

Same as the corresponding `new_V_*germline.tsv` **plus** filter annotations (all candidates, including filtered-out rows).

| Column | Type | Meaning |
|--------|------|---------|
| *(all candidates columns)* | | Same as in `candidates.tsv` |
| `is_filtered` | int | Number of filter reasons that applied (0 = kept) |
| `why_filtered` | string | Semicolon-separated reasons (e.g. `too_low_CDR3s_exact`, `xmap_ratio=...,other=...`) |
| `whitelist_diff` | int | 0 if consensus is in whitelist, -1 otherwise |

Filtered rows have `is_filtered > 0` and non-empty `why_filtered`.

---

## 8. `clonotypes.tsv`

One row per **clonotype**: representative row from the filtered table plus clonotype metadata. Columns are those of the **filtered table** (e.g. `v_call`, `j_call`, `cdr3`, `V_SHM`, …) **plus**:

| Column | Type | Meaning |
|--------|------|---------|
| `clonotype_id` | int | 1-based clonotype index |
| `clonotype_size` | int | Number of sequences (rows) in this clonotype |
| `CDR3_nt_mindiffrate` | float | When V_SHM is available: for the cluster, take the member with lowest V_SHM as reference; then 100 × (edit distance from representative CDR3 to reference CDR3) / length(reference CDR3). NaN if not computed (e.g. singleton or ref above threshold). |

The representative row is the one with **lowest V_SHM** in the clonotype (least mutated). Rows are a subset of the filtered table, so all other columns (v_call, j_call, cdr3, etc.) are unchanged from that table.

---

## 9. J discovery table (internal)

[`discover_j_genes`](@ref) returns a DataFrame with columns:

| Column | Type | Meaning |
|--------|------|---------|
| `name` | string | J allele name (database or unique name with sequence suffix) |
| `source` | string | J gene name |
| `exact` | int | Number of exact (zero-error) J alignments supporting this candidate |
| `CDR3s_exact` | int | Unique CDR3s among those exact matches |
| `consensus` | string | Consensus J sequence |

This table is not written as a standalone TSV in the pipeline; it is used to produce `new_J.fasta`.

---

Formulas above use the same logic as the code (see the Augment module for error/coverage/position helpers; V_SHM = 100 − identity). For AIRR base columns, see the [AIRR Rearrangement specification](https://docs.airr-community.org/).
