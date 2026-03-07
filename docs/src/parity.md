# Parity Testing

IgDiscover.jl includes Docker-based infrastructure that runs both the Python and Julia implementations on the same data and compares outputs.

## Running Parity Tests

```bash
make build    # Build Docker image
make test     # Run with synthetic reads (200 reads)
make test LIMIT=1000  # More reads
make test-reads READS=/path/to/reads.fasta.gz  # Real data
```

## What is Compared

- `filtered.tsv.gz` — Column-level comparison of preprocessed assignments
- `new_V_germline.fasta` — Discovered V gene sequences
- `candidates.tab` — Discovery statistics

## CI Integration

Parity tests run automatically on push to `main` via GitHub Actions.
