# IgDiscover.jl

A Julia port of [IgDiscover](https://gitlab.com/gkhlab/igdiscover22) for individualized V gene database construction from high-throughput sequencing data.

## Overview

IgDiscover.jl analyzes antibody repertoires and discovers new V genes from high-throughput sequencing reads. It supports heavy chains (VH), kappa (VK), and lambda (VL) light chains.

The package produces **identical results** to the Python implementation, verified through automated parity testing.

## Quick Start

```julia
using IgDiscover

# Initialize analysis directory
init_analysis("my_analysis", "path/to/database", "reads.fasta.gz")

# Run the pipeline
run_pipeline("my_analysis")
```

## Documentation

- **[How the algorithm works](@ref)** — Step-by-step: preprocessing, IgBLAST, augment, filter, discovery, germline filter, iteration, final run.
- **Pipeline** — File layout and iteration flow.
- **[Output formats](@ref)** — All output files and column definitions (TSV/JSON), with meanings and formulas.
- **[Configuration](@ref)** — TOML options and defaults.
- **Modules** — API by component (IgBLAST, grouping, augment, DNA, I/O, CDR3, alignment, clustering, discovery, J discovery, germline filter, clonotypes).
- **[API Reference](@ref)** — Exported functions and types.

## Features

- **Germline V gene discovery** — iterative consensus-based discovery with clustering
- **J gene discovery** — identifies novel J alleles
- **Preprocessing** — AIRR-format augmentation, filtering, IMGT database sanitization
- **Clonotype calling** — V+J+CDR3 single-linkage clustering
- **PCR bias correction** — barcode-based grouping with consensus
- **Configurable** — TOML-based configuration with sensible defaults
- **Fast** — threaded IgBLAST execution, pre-allocated edit distance buffers, PrecompileTools integration

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/IgDiscover.jl")
```

### External dependencies

IgDiscover.jl requires the same external tools as the Python version:

| Tool | Purpose | Install via |
|------|---------|-------------|
| `igblastn` | V/D/J gene assignment | conda: `bioconda::igblast` |
| `makeblastdb` | BLAST database creation | (bundled with igblast) |
| `muscle` | Multiple sequence alignment | conda: `bioconda::muscle` |

## Citation

> Corcoran et al. *Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity.* Nature Communications 7:13642 (2016)
