# Pipeline

The pipeline runs iteratively: each iteration discovers V gene candidates, filters them, and uses the result as the database for the next iteration.

## Steps

1. **IgBLAST** — Assign V/D/J genes to reads
2. **Augment** — Add IgDiscover-specific columns (errors, coverage, CDR3)
3. **Filter** — Apply preprocessing thresholds
4. **Discover** — Find V gene candidates via consensus
5. **Germline Filter** — Apply per-entry and pairwise filters
6. **J Discovery** — Discover J gene alleles (iteration 1 only)

## API

```@docs
init_analysis
run_pipeline
```
