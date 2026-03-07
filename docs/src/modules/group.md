# Read grouping

PCR bias correction: group reads by UMI (barcode) and optionally by pseudo-CDR3 or real CDR3, then output one representative sequence per group (either the first read or a consensus when the group is large enough).

```@autodocs
Modules = [IgDiscover]
Pages = ["group.jl"]
```
