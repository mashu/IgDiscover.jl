# IgDiscover.jl

Julia port of [IgDiscover](https://gitlab.com/gkhlab/igdiscover22), a tool for analyzing antibody repertoires and discovering new V genes from high-throughput sequencing reads.

## Design Principles

- **Result parity** with Python igdiscover, particularly `new_V_germline.tsv` and `filtered.tsv.gz`
- **Same external tools**: IgBLAST, MUSCLE, PEAR — no reimplementation of bioinformatics tools
- **Julia idioms**: multiple dispatch, concrete parametric types, functors, TOML config
- **No anti-patterns**: zero `try/catch`, `Any`, `isa()`, `typeof()`, `_`-prefixed functions
- **Performance**: thread-local buffers for hot-path edit distance, precomputed hash indices for exact-match lookups, zero-allocation nucleotide counting in consensus

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/IgDiscover.jl")
```

### External Dependencies

These must be in your `PATH`:

| Tool | Purpose | Install |
|------|---------|---------|
| `igblastn` | V/D/J gene assignment | [NCBI IgBLAST](https://ncbi.github.io/igblast/) |
| `makeblastdb` | Build BLAST databases | Ships with IgBLAST |
| `muscle` | Multiple sequence alignment | [MUSCLE](https://drive5.com/muscle/) |
| `pear` | Paired-end read merging | [PEAR](https://github.com/tseemann/PEAR) |

Consider using [IgBLAST.jl](https://github.com/mashu/IgBLAST.jl) for automatic IgBLAST binary management.

## Quick Start

### Initialize an Analysis

```julia
using IgDiscover

init_analysis("my_analysis",           # output directory
              "path/to/database",      # directory with V.fasta, D.fasta, J.fasta
              "path/to/reads.fasta.gz") # merged reads
```

### Edit Configuration

```toml
# my_analysis/igdiscover.toml
iterations = 1
seed = 1

[germline_filter]
unique_cdr3s = 5
unique_js = 3
cluster_size = 50
```

### Run the Pipeline

```julia
run_pipeline("my_analysis")
```

**Output files** (per iteration, e.g. `iteration-01/`):

| File | Description |
|------|-------------|
| `airr.tsv.gz` | IgBLAST AIRR-format assignments |
| `assigned.tsv.gz` | Augmented table with IgDiscover columns |
| `filtered.tsv.gz` | Rows passing preprocessing filter |
| `candidates.tab` | V gene candidates from discovery |
| `new_V_pregermline.fasta` | Pre-germline filtered sequences |
| `new_V_germline.fasta` | **Final discovered V germline sequences** |
| `new_J.fasta` | Discovered J genes |

### IMGT Database Handling

IgDiscover.jl automatically sanitizes IMGT-format databases during initialization:
- Extracts allele names from pipe-delimited headers (e.g. `M99641|IGHV1-18*01|...` → `IGHV1-18*01`)
- Removes IMGT gap dots from sequences (e.g. `cagg.ttcag...tctgg` → `CAGGTTCAGTCTGG`)

This replaces the need for `edit_imgt.pl` or manual preprocessing.

## Architecture

```
src/
├── IgDiscover.jl     # Module definition
├── dna.jl            # DNA utilities (translate, edit_distance with thread-local buffers)
├── config.jl         # TOML configuration with concrete types
├── io.jl             # FASTA/TSV I/O, IMGT sanitization
├── cdr3.jl           # CDR3 detection (ported from species.py)
├── alignment.jl      # Multiple alignment, consensus (NucleotideCounter), affine alignment
├── clustering.jl     # Hierarchical + single-linkage clustering
├── group.jl          # PCR bias correction (UMI/barcode grouping)
├── igblast.jl        # IgBLAST wrapper (calls external igblastn)
├── augment.jl        # AIRR table augmentation with IgDiscover columns
├── tablefilter.jl    # Preprocessing filter
├── discovery.jl      # V gene candidate discovery (core algorithm)
├── jdiscovery.jl     # J gene candidate discovery
├── germlinefilter.jl # Germline/pre-germline filter (dispatch-based)
├── rename.jl         # Gene renaming (canonical names from database)
└── pipeline.jl       # Full pipeline orchestration
```

### Key Design Decisions

**Germline filters use abstract dispatch:**
```julia
abstract type CandidateFilter end
struct CrossMappingFilter <: CandidateFilter ... end
struct ExactRatioFilter <: CandidateFilter ... end

should_discard(f::CrossMappingFilter, ref, cand, same_gene) = ...
should_discard(f::ExactRatioFilter, ref, cand, same_gene) = ...
```

**Thread-local edit distance buffers (zero allocation in hot loops):**
```julia
# Pre-allocated per-thread buffers avoid GC pressure during O(n²) pairwise comparisons
edit_distance("ACGT", "AGGT"; maxdiff=1)  # 0 allocations after warmup
```

**Zero-allocation consensus counting:**
```julia
# NucleotideCounter replaces Dict{Char,Int} in the consensus inner loop
# Fixed-size struct with fields for A, C, G, T, N, gap
```

**Precomputed hash indices in discovery:**
```julia
# O(1) exact-match lookups instead of O(n) linear scans
v_no_cdr3_index = Dict{String,Vector{Int}}()  # sequence → row indices
```

## Testing

### Unit Tests

```julia
using Pkg
Pkg.test("IgDiscover")
```

### Docker Parity Test

The project includes a Dockerfile that installs both Python igdiscover (via bioconda)
and Julia IgDiscover.jl, runs both on the same data, and compares outputs.

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

The Docker image uses `mambaorg/micromamba` and installs:
- Python igdiscover 0.15.1 (bioconda: includes igblast 1.17, muscle 3.8, pear 0.9.6)
- Julia 1.11.3 LTS with IgDiscover.jl precompiled

### Manual Parity Comparison

```bash
julia --project=. test/run_parity_test.jl /path/to/python/analysis /path/to/julia/analysis
```

## Roadmap

- [x] Core germline discovery pipeline
- [x] Preprocessing filter (filtered.tsv.gz)
- [x] Germline filter (new_V_germline.fasta)
- [x] CDR3 detection from V/J reference
- [x] TOML configuration
- [x] Thread-local edit distance buffers
- [x] Precomputed hash indices for discovery
- [x] J gene discovery
- [x] Gene renaming
- [x] PCR bias correction (barcode grouping)
- [x] IMGT database sanitization (replaces edit_imgt.pl)
- [x] Docker parity test infrastructure
- [x] Synthetic test data generation
- [ ] Clonotypes subcommand
- [ ] Clonoquery subcommand
- [ ] Plot alleles
- [ ] SnoopCompile + PackageCompiler binary bundle

## Citation

> Corcoran et al. *Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity.* Nature Communications 7:13642 (2016)

## License

MIT
