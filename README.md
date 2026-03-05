# IgDiscover.jl

Julia port of [IgDiscover](https://gitlab.com/gkhlab/igdiscover22), a tool for analyzing antibody repertoires and discovering new V genes from high-throughput sequencing reads.

## Design Principles

- **Result parity** with Python igdiscover, particularly `new_V_germline.tsv` and `filtered.tsv.gz`
- **Same external tools**: IgBLAST, MUSCLE, PEAR — no reimplementation of bioinformatics tools
- **Julia idioms**: multiple dispatch, concrete parametric types, functors, TOML config
- **No anti-patterns**: zero `try/catch`, `Any`, `isa()`, `typeof()`, `_`-prefixed functions

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

This creates:
```
my_analysis/
├── igdiscover.toml     # configuration (edit before running)
├── database/
│   ├── V.fasta
│   ├── D.fasta
│   └── J.fasta
└── reads.fasta.gz
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

### Step-by-Step Usage

```julia
using IgDiscover

# Load config
config = load_config("igdiscover.toml")

# Run IgBLAST
run_igblast_on_fasta("database/", "reads.fasta.gz", "airr.tsv.gz";
    species="human", sequence_type="Ig")

# Augment with IgDiscover columns
df = read_assignments("airr.tsv.gz")
augmented = augment_table(df, "database/")
write_table_gzipped("assigned.tsv.gz", augmented)

# Filter
filtered, stats = filter_table(augmented, config.preprocessing_filter)
write_table_gzipped("filtered.tsv.gz", filtered)

# Discover candidates
database = read_fasta_dict("database/V.fasta")
candidates = discover_germline(filtered, database, config)
write_table("candidates.tab", candidates)

# Germline filter
filtered_candidates, annotated = germline_filter!(candidates, config.germline_filter)
germline_filter_to_fasta(filtered_candidates, "new_V_germline.fasta")
```

## Architecture

```
src/
├── IgDiscover.jl     # Module definition
├── dna.jl            # DNA utilities (translate, edit_distance, etc.)
├── config.jl         # TOML configuration with concrete types
├── io.jl             # FASTA/TSV I/O
├── cdr3.jl           # CDR3 detection (ported from species.py)
├── alignment.jl      # Multiple alignment, consensus, affine alignment
├── clustering.jl     # Hierarchical + single-linkage clustering
├── igblast.jl        # IgBLAST wrapper (calls external igblastn)
├── augment.jl        # AIRR table augmentation with IgDiscover columns
├── tablefilter.jl    # Preprocessing filter (V/J coverage, evalue)
├── discovery.jl      # V gene candidate discovery (core algorithm)
├── germlinefilter.jl # Germline/pre-germline filter (dispatch-based)
└── pipeline.jl       # Full pipeline orchestration
```

### Key Design Decisions

**Germline filters use abstract dispatch:**
```julia
abstract type CandidateFilter end
struct CrossMappingFilter <: CandidateFilter ... end
struct ExactRatioFilter <: CandidateFilter ... end

# Each filter dispatches on its concrete type:
should_discard(f::CrossMappingFilter, ref, cand, same_gene) = ...
should_discard(f::ExactRatioFilter, ref, cand, same_gene) = ...
```

**Configuration uses immutable concrete structs:**
```julia
struct GermlineFilterCriteria
    unique_cdr3s::Int
    unique_js::Int
    # ... all fields are concrete types
end
```

## Testing

```julia
using Pkg
Pkg.test("IgDiscover")
```

### Parity Testing

To verify exact parity with Python igdiscover:

1. Run Python igdiscover on test data (e.g., ERR1760498)
2. Run Julia IgDiscover.jl on the same data
3. Compare outputs:

```bash
julia --project=. scripts/run_parity_test.jl /path/to/python/analysis /path/to/julia/analysis
```

## Roadmap

- [x] Core germline discovery pipeline
- [x] Preprocessing filter (filtered.tsv.gz)
- [x] Germline filter (new_V_germline.fasta)
- [x] CDR3 detection from V/J reference
- [x] TOML configuration
- [x] Parity test infrastructure
- [ ] J gene discovery (discoverjd)
- [ ] Clonotypes subcommand
- [ ] Clonoquery subcommand
- [ ] Plot alleles
- [ ] SnoopCompile + PackageCompiler binary bundle

## Citation

If you use IgDiscover, please cite:

> Corcoran et al. *Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity.* Nature Communications 7:13642 (2016)

## License

MIT
