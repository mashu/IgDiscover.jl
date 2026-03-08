# IgDiscover.jl

[![CI](https://github.com/mashu/IgDiscover.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mashu/IgDiscover.jl/actions/workflows/CI.yml)
[![Release](https://github.com/mashu/IgDiscover.jl/actions/workflows/Release.yml/badge.svg)](https://github.com/mashu/IgDiscover.jl/actions/workflows/Release.yml)
[![codecov](https://codecov.io/gh/mashu/IgDiscover.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/IgDiscover.jl)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/IgDiscover.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://mashu.github.io/IgDiscover.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Parity](https://img.shields.io/badge/parity-verified-brightgreen.svg)](#parity-testing)

A reimplementation in Julia of [IgDiscover](https://gitlab.com/gkhlab/igdiscover22) for individualized V gene database construction from high-throughput antibody sequencing data.

IgDiscover.jl produces **identical results** to the Python implementation, verified through automated Docker-based parity testing. It is designed as an idiomatic Julia system using parametric types, multiple dispatch, callable structs, and compile-time type resolution — not a mechanical translation of the Python codebase.

## Quick Start

```julia
using IgDiscover

# Initialize analysis directory
init_analysis("my_analysis", "path/to/database", "reads.fasta.gz")

# Run the full pipeline
run_pipeline("my_analysis")
```

Or from the command line:

```bash
igdiscover init my_analysis database/ reads.fasta.gz
igdiscover run my_analysis
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

### Core Pipeline

The main germline discovery pipeline runs iteratively: IgBLAST assignment → augmentation → filtering → candidate discovery → germline filtering → gene renaming, with J discovery in the first iteration.

- **Germline V gene discovery** — iterative consensus with SHM-window and cluster-based candidate generation
- **J gene discovery** — novel J allele identification with cross-mapping filter
- **Preprocessing** — AIRR augmentation, coverage/evalue filtering, IMGT database sanitization
- **PCR bias correction** — barcode grouping with consensus
- **Gene renaming** — canonical allele naming based on database proximity
- **TOML configuration** — sensible defaults, easy to override

### Analysis Subcommands

All of the original IgDiscover subcommands have been ported:

| Subcommand | Description |
|------------|-------------|
| `clonotypes` | V+J+CDR3 single-linkage clonotype clustering |
| `clonoquery` | Query a reference table by clonotype similarity |
| `count` | Gene expression counting with allele ratio filtering |
| `plotalleles` | Allele co-expression matrix computation |
| `haplotype` | Haplotype phasing from allele co-occurrences |
| `upstream` | UTR/leader consensus per V gene |
| `clusterplot` | Per-gene sequence clustering with distance heatmaps |
| `dbdiff` | Compare two FASTA gene databases |
| `dereplicate` | Sequence deduplication with barcode removal |
| `merge` | Paired-end read merging (PEAR wrapper) |
| `union` | Compute union of sequences across FASTA files |

### Performance

- Threaded IgBLAST execution (chunked, one process per chunk)
- Threaded O(n²) pairwise distance matrix computation
- Thread-local pre-allocated edit distance buffers (zero-allocation hot loops)
- `NucleotideCounter` struct replacing `Dict{Char,Int}` in consensus
- MUSCLE v3/v5 auto-detection with version-appropriate flags
- PrecompileTools workload for reduced first-call latency
- Optional PackageCompiler sysimage for near-instant startup

### Architecture

The codebase is organized around Julia's type system and multiple dispatch:

- **Callable structs (functors)** for pipeline steps: `ClonotypeCaller`, `ClonoQuery`, `AlleleUsage`, `HaplotypeAnalyzer`, `ClusterPlotter`, `ExpressionCounter`, `DatabaseComparator`, `UpstreamAnalyzer`, `Dereplicator`, `ReadMerger`, `SequenceUnion`
- **Abstract type hierarchies** with dispatch: `CandidateFilter` → `IdenticalSequenceFilter`, `CrossMappingFilter`, `ClonotypeRatioFilter`, etc.
- **Parametric types**: `UnionFind{T}`, `EditDistanceBuffers`, `MuscleAligner{V}`, `ReadMerger{P}`
- **No runtime type checks**: no `Any`, `isa`, `typeof`, `eltype`, or try-catch flow control
- **Clean dispatch**: `apply_filters(table, ::AlleleUsage)`, `apply_filters(table, ::HaplotypeAnalyzer)`, `apply_filters(table, ::ExpressionCounter)` — same function name, different behavior via type

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

IgDiscover.jl produces results identical to the Python implementation, verified through Docker-based automated testing. Both `filtered.tsv.gz` and `new_V_germline.fasta` are compared for bit-identity.

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

The Docker image installs both Python igdiscover 0.15.1 (via bioconda, including igblast, muscle, pear) and Julia IgDiscover.jl, runs both on identical data, and compares outputs.

## Binaries (no Julia required)

Pre-built standalone binaries are available on [GitHub Releases](https://github.com/mashu/IgDiscover.jl/releases) for Linux (x86_64), macOS (ARM64), and Windows (x86_64). They are built automatically on every tagged release.

1. Download the archive for your platform from the [latest release](https://github.com/mashu/IgDiscover.jl/releases/latest).
2. Extract it and run `bin/igdiscover` (or `bin/igdiscover.exe` on Windows).
3. Ensure [external dependencies](#external-dependencies) (IgBLAST, MUSCLE, optionally PEAR) are installed.

```bash
./bin/igdiscover --help
./bin/igdiscover init my_analysis database/ reads.fasta.gz
./bin/igdiscover run my_analysis
```

## Building a System Image

For near-instant startup (useful in production or interactive use):

```bash
# Install build dependencies (one-time; PackageCompiler is in [extras])
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(; target="build")'

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
├── utils.jl           # Shared utilities: tallies, gene_family, ensure_column!, from_toml
├── dna.jl             # DNA utilities: translate, edit_distance, hamming_distance
├── config.jl          # TOML configuration types and loading
├── io.jl              # FASTA/FASTQ I/O, IMGT sanitization, TSV I/O
├── cdr3.jl            # CDR3 start/end detection from V/J reference
├── alignment.jl       # MUSCLE wrapper (v3/v5), consensus, affine alignment
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
├── clonoquery.jl      # Clonoquery — query reference table by clonotype
├── plotalleles.jl     # Allele co-expression matrix
├── upstream.jl        # UTR/leader consensus per V gene
├── count.jl           # Gene expression counting
├── dbdiff.jl          # FASTA database comparison
├── haplotype.jl       # Haplotype phasing from allele co-occurrences
├── clusterplot.jl     # Per-gene sequence clustering + distance heatmap
├── dereplicate.jl     # Sequence deduplication with barcode removal
├── merge.jl           # Paired-end read merging (PEAR/FLASH wrapper)
├── union.jl           # Compute union of FASTA sequences (prefix merging)
├── pipeline.jl        # Pipeline orchestration with file-level caching
└── cli.jl             # CLI (ArgParse) + PackageCompiler entry point
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
- [x] Clonoquery subcommand
- [x] Expression counting (`count` subcommand)
- [x] Allele co-expression matrix (`plotalleles` subcommand)
- [x] Haplotype analysis (`haplotype` subcommand)
- [x] UTR/leader consensus (`upstream` subcommand)
- [x] Database comparison (`dbdiff` subcommand)
- [x] Cluster visualization (`clusterplot` subcommand)
- [x] Sequence dereplication (`dereplicate` subcommand)
- [x] Paired-end read merging (`merge` subcommand)
- [x] Sequence union across FASTA files (`union` subcommand)
- [x] Per-iteration expression counts and exact tables
- [x] Per-iteration clusterplot data export
- [x] GitHub CI + codecov
- [x] Documenter.jl documentation
- [x] PackageCompiler build scripts
- [x] MUSCLE v3/v5 auto-detection
- [ ] PackageCompiler sysimage in CI artifacts
- [ ] Allele ratio plots (data export; plotting deferred to downstream tools)
- [ ] Error histograms (data export; plotting deferred to downstream tools)
- [ ] Dendrograms (data export; plotting deferred to downstream tools)

## Citation

> Corcoran et al. *Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity.* Nature Communications 7:13642 (2016)

## License

MIT
