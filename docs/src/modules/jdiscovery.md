# J Discovery

J gene discovery: find novel J alleles from the filtered assignment table (clustering, consensus, allele/cross-mapping filters), then write `new_J.fasta`. Used in the first iteration; results can be propagated to later iterations and the final run.

```@docs
IgDiscover.discover_j_genes
IgDiscover.discover_j_to_fasta
IgDiscover.filter_j_alleles
IgDiscover.filter_j_cross_mapping
```
