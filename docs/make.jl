using Documenter
using IgDiscover

DocMeta.setdocmeta!(IgDiscover, :DocTestSetup, :(using IgDiscover); recursive=true)

makedocs(;
    modules=[IgDiscover],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com>",
    sitename="IgDiscover.jl",
    format=Documenter.HTML(;
        canonical="https://mashu.github.io/IgDiscover.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Algorithm" => "algorithm.md",
        "Output formats" => "output_formats.md",
        "Configuration" => "config.md",
        "Pipeline" => "pipeline.md",
        "Modules" => [
            "IgBLAST" => "modules/igblast.md",
            "Read grouping" => "modules/group.md",
            "Augment" => "modules/augment.md",
            "DNA Utilities" => "modules/dna.md",
            "FASTA I/O" => "modules/io.md",
            "CDR3 Detection" => "modules/cdr3.md",
            "Alignment" => "modules/alignment.md",
            "Clustering" => "modules/clustering.md",
            "Discovery" => "modules/discovery.md",
            "J Discovery" => "modules/jdiscovery.md",
            "Germline Filter" => "modules/germlinefilter.md",
            "Clonotypes" => "modules/clonotypes.md",
        ],
        "Parity Testing" => "parity.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/mashu/IgDiscover.jl",
    devbranch="main",
)
