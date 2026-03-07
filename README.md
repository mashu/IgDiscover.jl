# IgDiscover.jl

[![CI](https://github.com/mashu/IgDiscover.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mashu/IgDiscover.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/mashu/IgDiscover.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgDiscover.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/IgDiscover.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgDiscover.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Parity](https://img.shields.io/badge/parity-verified-brightgreen.svg)](#parity-testing)

A Julia port of [IgDiscover](https://gitlab.com/gkhlab/igdiscover22) for individualized V gene database construction from high-throughput antibody sequencing data.

## Quick Start

```julia
using IgDiscover

# Initialize analysis directory
init_analysis("my_analysis", "path/to/database", "reads.fasta.gz")

# Run the full pipeline
run_pipeline("my_analysis")
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/IgDiscover.jl")
```

### External Dependencies

IgDiscover.jl relies on the same external tools as the Python version:

| Tool | Purpose | Install |
|------|---------|---------|
| `igblastn` + `makeblastdb` | V/D/J gene assignment | `conda install -c bioconda igblast` |
| `muscle` | Multiple sequence alignment | `conda install -c bioconda muscle` |
| `pear` (optional) | Read merging | `conda install -c bioconda pear` |

## Features

- **Germline V gene discovery** — iterative consensus with SHM-window and cluster-based candidate generation
- **J gene discovery** — novel J allele identification with cross-mapping filter
- **Preprocessing** — AIRR augmentation, coverage/evalue filtering, IMGT database sanitization
- **Clonotype calling** — V+J+CDR3 single-linkage clustering
- **PCR bias correction** — barcode grouping with consensus
- **Gene renaming** — canonical allele naming based on database proximity
- **TOML configuration** — with sensible defaults, easy to override

### Performance

- Threaded IgBLAST execution (chunked, one process per chunk)
- Threaded O(n²) pairwise distance matrix computation
- Thread-local pre-allocated edit distance buffers (zero-allocation hot loops)
- `NucleotideCounter` struct replacing `Dict{Char,Int}` in consensus
- PrecompileTools workload for reduced first-call latency
- Optional PackageCompiler sysimage for near-instant startup

## Configuration

The pipeline uses TOML configuration. See [`config/defaults.toml`](config/defaults.toml) for all options:

```toml
species = ""
sequence_type = "Ig"
iterations = 1
limit = 0                     # 0 = process all reads
consensus_threshold = 60.0

[preprocessing_filter]
v_coverage = 90.0
j_coverage = 60.0
v_evalue = 1e-3

[germline_filter]
unique_cdr3s = 5
unique_js = 3
cluster_size = 50
cross_mapping_ratio = 0.02
```

## Parity Testing

IgDiscover.jl produces results identical to the Python implementation, verified through Docker-based automated testing.

```bash
# Build the parity test image
make build

# Run with synthetic reads (fast, ~200 reads)
make test

# Run with real reads
make test-reads READS=/path/to/reads.fasta.gz LIMIT=1000

# Interactive debugging
make test-interactive
```

The Docker image installs both Python igdiscover 0.15.1 (via bioconda, including igblast, muscle, pear) and Julia IgDiscover.jl, runs both on identical data, and compares `filtered.tsv.gz` and `new_V_germline.fasta` outputs.

## Building a System Image

For near-instant startup (useful in production or interactive use):

```bash
# Install PackageCompiler (one-time)
julia -e 'using Pkg; Pkg.add("PackageCompiler")'

# Build system image (~5 min)
julia --project=. build/build_sysimage.jl

# Use it
julia --sysimage=build/IgDiscover.so --project=. bin/igdiscover.jl run

# Or build a standalone application
julia --project=. build/build_sysimage.jl --app
```

## Project Structure

```
src/
├── IgDiscover.jl      # Module definition + PrecompileTools workload
├── dna.jl             # DNA utilities: translate, edit_distance, hamming_distance
├── config.jl          # TOML configuration types and loading
├── io.jl              # FASTA/FASTQ I/O, IMGT sanitization, TSV I/O
├── cdr3.jl            # CDR3 start/end detection from V/J reference
├── alignment.jl       # MUSCLE wrapper, consensus, affine alignment
├── clustering.jl      # Hierarchical + single-linkage clustering, UnionFind
├── group.jl           # PCR bias correction (barcode grouping)
├── igblast.jl         # IgBLAST wrapper (threaded chunked execution)
├── augment.jl         # AIRR table augmentation (errors, coverage, CDR3)
├── tablefilter.jl     # Preprocessing filter (→ filtered.tsv.gz)
├── discovery.jl       # V gene candidate discovery
├── jdiscovery.jl      # J gene discovery
├── germlinefilter.jl  # Germline filter (abstract CandidateFilter + dispatch)
├── rename.jl          # Gene renaming to canonical names
├── clonotypes.jl      # Clonotype calling (V+J+CDR3 clustering)
└── pipeline.jl        # Pipeline orchestration
```

## Roadmap

- [x] Core germline discovery pipeline
- [x] Preprocessing filter (`filtered.tsv.gz`)
- [x] Germline filter (`new_V_germline.fasta`)
- [x] CDR3 detection from V/J reference
- [x] TOML configuration
- [x] Thread-local edit distance buffers
- [x] Precomputed hash indices for discovery
- [x] J gene discovery
- [x] Gene renaming
- [x] PCR bias correction (barcode grouping)
- [x] IMGT database sanitization
- [x] Docker parity test infrastructure
- [x] Synthetic test data generation
- [x] PrecompileTools startup workload
- [x] Clonotype calling
- [x] GitHub CI + codecov
- [x] Documenter.jl documentation
- [x] PackageCompiler build scripts
- [ ] Clonoquery subcommand
- [ ] Plot alleles (data export for plotting)
- [ ] Allele ratio plots
- [ ] Upstream/downstream analysis
- [ ] PackageCompiler sysimage in CI artifacts

## Citation

> Corcoran et al. *Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity.* Nature Communications 7:13642 (2016)

## License

MIT
