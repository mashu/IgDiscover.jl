# How the algorithm works

This page explains IgDiscover’s workflow in a structured, step-by-step way: what each stage does and why it matters for building an individualized V gene database.

---

## Goal

IgDiscover builds **individualized germline V (and J) gene databases** from high-throughput sequencing of antibody or TCR repertoires. Public references (e.g. IMGT) may miss alleles or genes present in the donor. By assigning reads to genes, clustering by V assignment, and building consensus sequences, IgDiscover discovers **novel V (and J) alleles** and produces a donor-specific database suitable for downstream assignment and repertoire analysis.

---

## Pipeline at a glance

The pipeline has two layers:

1. **Optional read-level preprocessing** (before iteration): barcode-based grouping and/or RACE-G trimming.
2. **Iterative discovery**: for each iteration, (a) assign reads to the current V/D/J database, (b) augment and filter the assignment table, (c) discover V candidates from clusters, (d) apply germline filters, (e) optionally discover J genes (first iteration only). The next iteration uses the new V (and optionally J) database so that previously unmapped reads can be assigned and new alleles can be refined.

After the last iteration, a **final** run reassigns reads to the final database and produces the output tables and clonotype calls.

---

## Step 1: Read preprocessing (optional)

**What it does**

- **Barcode grouping** (when `barcode_length != 0`): Groups reads by UMI (barcode). Within each barcode, reads can be sub-grouped by *pseudo-CDR3* or *real CDR3* (if available). For each group, one representative is written: either the first read or a consensus sequence when the group is large enough and the consensus is unambiguous. This reduces PCR duplicates and avoids inflating support for a single molecule.
- **RACE-G trimming** (when `race_g` is true and no barcode): Removes leading G nucleotides (RACE protocol artifact) from each read.

**Why it matters**

- Grouping reduces PCR bias and keeps the assignment step focused on unique molecules.
- Clean 5′ ends improve V alignment and assignment.

Relevant functions: [`group_reads`](@ref), [`IgDiscover.trim_reads_race_g`](@ref). See the Read grouping module.

---

## Step 2: IgBLAST assignment

**What it does**

- Builds BLAST databases from the current V, D, and J FASTA files in the iteration (or final) directory.
- Runs `igblastn` on the read FASTA (or grouped/trimmed FASTA) and writes AIRR-format tabular output (e.g. `airr.tsv.gz`): one row per read with columns such as `v_call`, `j_call`, `d_call`, alignment coordinates, and sequence alignments.

**Why it matters**

- IgBLAST provides the initial V/D/J assignment and alignment that everything else builds on. The same reads are re-run in later iterations with an updated V (and optionally J) database so that new alleles get assigned and refined.

Relevant functions: [`run_igblast_on_fasta`](@ref), [`IgDiscover.make_vdj_blastdb`](@ref), [`make_blastdb`](@ref). See the IgBLAST module.

---

## Step 3: Augment

**What it does**

- Takes the raw AIRR table and adds IgDiscover-specific columns using the alignment strings and the reference databases:
  - **Error counts**: V/D/J alignment error counts (mismatches + gaps).
  - **Coverage**: Fraction of the reference (V/D/J) covered by the alignment.
  - **CDR3**: Nucleotide CDR3 and its boundaries (using locus-specific rules and alignment positions).
  - **Header-derived fields**: From FASTA headers (e.g. `size=`, `barcode=`) when present.

**Why it matters**

- Filtering and discovery rely on these columns (e.g. V coverage, no internal stop, CDR3 present). Augment is the place that turns raw IgBLAST output into the table used by the rest of the pipeline.

Relevant functions: [`augment_table`](@ref), [`parse_header`](@ref), [`IgDiscover.count_alignment_errors`](@ref), [`IgDiscover.alignment_coverage`](@ref), [`IgDiscover.query_position`](@ref). See the Augment module.

---

## Step 4: Preprocessing filter

**What it does**

- Keeps only rows that have both V and J assignment, no internal stop codon, and pass configurable thresholds (e.g. V/J coverage, V E-value). Counts at each stage are recorded as [`FilterStats`](@ref).

**Why it matters**

- Low-quality or incomplete assignments would add noise to clustering and consensus. This step keeps the table focused on confident, full-length V–J assignments.

Relevant functions: [`filter_table`](@ref). Configuration: [`PreprocessingFilter`](@ref).

---

## Step 5: V gene discovery

**What it does**

- Groups the filtered table by **V gene assignment** (`v_call`).
- For each V gene with enough reads:
  - Partitions reads into **clusters** (e.g. by sequence similarity / windows).
  - Builds a **consensus** sequence per cluster (multiple sequence alignment + majority-rule consensus).
  - For each consensus, computes summary stats: cluster size, unique J/CDR3/D, clonotype count, CDR3 ratios, comparison to the current database sequence (edit distance, HGVS-like changes), etc. These become **V gene candidates** (see the Discovery module; written as `candidates.tsv`).
- Candidates that exactly match the current database sequence are re-included as “expressed” if they have enough support.

**Why it matters**

- Consensus turns many reads into one representative sequence per cluster, which is the basis for new alleles. The stats (e.g. unique CDR3s, clonotypes) are used later by the germline filter to decide which candidates to keep.

Relevant functions: [`discover_germline`](@ref). See the Discovery and Alignment modules (consensus, multialign).

---

## Step 6: Germline filter

**What it does**

- Two stages (both configurable): **pre-germline** and **germline**.
- **Per-entry**: Drops candidates that fail criteria such as minimum unique CDR3s, minimum unique Js, minimum cluster size, maximum N count, etc.
- **Pairwise**: For pairs of similar candidates (e.g. same gene, small edit distance), applies filters (e.g. identical sequence, cross-mapping ratio, clonotype ratio, exact occurrence ratio, unique D ratio). The “weaker” candidate (e.g. lower support) can be discarded with a reason.
- **Whitelist**: Known sequences (e.g. from the starting database or a `whitelist.fasta`) can be marked so they are not discarded by these rules.
- Outputs: filtered candidate table, annotated table (with discard reasons), and a FASTA of accepted germline V sequences (`new_V_pregermline.fasta`, `new_V_germline.fasta`). Optionally **deduplicate by consensus** so that one row per unique sequence is written (matching Python behavior).

**Why it matters**

- Not every consensus is a true germline allele (e.g. PCR/sequencing errors, rare alleles with little support). The germline filter uses support, diversity (CDR3/J), and pairwise comparisons to keep likely true alleles and drop artifacts or low-support duplicates.

Relevant functions: [`germline_filter!`](@ref), [`deduplicate_by_consensus`](@ref), [`germline_filter_to_fasta`](@ref). Configuration: [`GermlineFilterCriteria`](@ref). See the Germline Filter module.

---

## Step 7: J gene discovery (first iteration only)

**What it does**

- Uses the filtered assignment table and current J database to find J gene candidates (novel alleles) in a similar spirit to V discovery: clustering, consensus, and filtering. Results are written as `new_J.fasta` (and related files). If configured, this updated J database can be **propagated** to later iterations and to the final run.

**Why it matters**

- Individualized J alleles improve assignment accuracy and CDR3 boundary detection, which in turn improves V discovery and downstream analyses.

Relevant functions: [`discover_j_genes`](@ref), [`discover_j_to_fasta`](@ref). Configuration: [`JDiscoveryConfig`](@ref). See the J Discovery module.

---

## Step 8: Iteration and final run

**What it does**

- **Iteration**: The next iteration uses the new V (and optionally J) database from the previous one. Reads are *re-run* through IgBLAST with this database, then augment → filter → discover → germline (and J discovery only in iteration 1). So each iteration can assign previously unmapped reads to newly discovered alleles and refine candidates.
- **Final**: After the last iteration, the pipeline copies the final V/D/J databases into a `final/` directory, re-runs IgBLAST (and augment/filter) with the final database unless results can be reused, then runs **clonotype calling** on the final filtered table: group sequences by V, J, and CDR3 similarity (e.g. Hamming distance ≤ threshold) and write a clonotype table.

**Why it matters**

- Iteration is what allows the database to grow and improve from the data. The final run produces the canonical assignment and clonotype outputs for the user.

Relevant functions: [`run_pipeline`](@ref), [`init_analysis`](@ref), [`call_clonotypes`](@ref). See the pipeline page and Clonotypes module.

---

## Summary diagram

```
Reads (FASTA)
    → [optional] group by barcode / trim RACE-G
    → IgBLAST (V/D/J assignment) → airr.tsv.gz
    → Augment (errors, coverage, CDR3) → assigned.tsv.gz
    → Preprocessing filter → filtered.tsv.gz
    → Discover V candidates (cluster + consensus) → candidates.tsv
    → Germline filter (per-entry + pairwise) → new_V_germline.fasta
    → [iteration 1 only] J discovery → new_J.fasta
    → [next iteration: repeat with new V/J DB] …
    → Final: assign with final DB, clonotype call → final/
```

For configuration options and file layout, see Configuration and the pipeline page.
