# IgBLAST

IgBLAST assigns V, D, and J gene segments to each read. This module wraps the external `igblastn` and `makeblastdb` programs: it builds BLAST databases from the germline FASTA files and runs IgBLAST on read chunks (optionally in parallel).

```@autodocs
Modules = [IgDiscover]
Pages = ["igblast.jl"]
```
