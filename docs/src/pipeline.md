# Pipeline

The pipeline runs **iteratively**: each iteration discovers V gene candidates, applies germline filters, and uses the resulting V (and optionally J) database for the next iteration. This section summarizes the flow and file layout; for a detailed explanation of each step, see [How the algorithm works](@ref).

## Iteration flow

For each iteration (e.g. `iteration-01`, `iteration-02`, …):

1. **IgBLAST** — Assign V/D/J genes to reads using the current database. Output: `iteration-NN/airr.tsv.gz`.
2. **Augment** — Add IgDiscover columns (errors, coverage, CDR3, header fields). Output: `iteration-NN/assigned.tsv.gz`.
3. **Filter** — Apply preprocessing thresholds (V+J assigned, no stop, coverage/evalue). Output: `iteration-NN/filtered.tsv.gz`, `iteration-NN/stats/filtered.json`.
4. **Discover** — Build V gene candidates from clusters and consensus. Output: `iteration-NN/candidates.tsv`.
5. **Germline filter** — Per-entry and pairwise filters; output `new_V_pregermline.fasta` and `new_V_germline.fasta` (and annotated `.tsv` files).
6. **J discovery** — Only in the first iteration: discover J alleles → `iteration-01/new_J.fasta`. Optional propagation to later iterations and final run.

The **next** iteration uses `new_V_pregermline.fasta` (or the previous iteration’s germline output) as the V database, so newly discovered alleles are used for re-assignment and refinement.

## Final run

After the last iteration, the pipeline:

- Copies the final V/D/J databases into `final/database/`.
- Re-runs IgBLAST (and augment/filter) with this database unless results can be reused.
- Runs **clonotype calling** on the final filtered table → `final/clonotypes.tsv`.

Outputs in `final/`: `airr.tsv.gz`, `assigned.tsv.gz`, `filtered.tsv.gz`, `database/V.fasta` (and D, J), `stats/`, and `clonotypes.tsv`. For column definitions and formulas, see [Output formats](@ref).

## Optional preprocessing (before iterations)

- If **barcode grouping** is enabled (`barcode_length != 0`), reads are grouped by UMI (and optionally by CDR3), and one representative per group is written (e.g. `grouped.fasta.gz`) before the first IgBLAST run.
- If **RACE-G trim** is enabled (`race_g` and no barcode), leading G nucleotides are trimmed (e.g. `trimmed.fasta.gz`).

## API

- [`init_analysis`](@ref) — Set up analysis directory with database, reads, and default config.
- [`run_pipeline`](@ref) — Run the full pipeline (preprocessing, iterations, final, clonotypes).

For all pipeline-related functions and types, see the [API Reference](@ref).
