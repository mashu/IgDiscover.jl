# Configuration

IgDiscover.jl uses TOML configuration files. A default configuration is written when you call `init_analysis`.

## Configuration File

```toml
species = ""
sequence_type = "Ig"        # "Ig" or "TCR"
iterations = 1
limit = 0                   # 0 = process all reads

[preprocessing_filter]
v_coverage = 90.0
j_coverage = 60.0
v_evalue = 1e-3

[germline_filter]
unique_cdr3s = 5
unique_js = 3
cluster_size = 50
# ... see defaults.toml for all options
```

## API

```@docs
Config
load_config
write_default_config
PreprocessingFilter
GermlineFilterCriteria
JDiscoveryConfig
```
