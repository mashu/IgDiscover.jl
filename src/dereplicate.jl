# Dereplicate — deduplicate sequences and remove barcodes
#
# Unlike the `group` module which computes consensus from barcode/CDR3 groups,
# this module simply removes duplicate sequences and strips barcodes.
# Barcodes can be at the 5' or 3' end (positive/negative barcode_length).

"""
    DereplicateParams(; barcode_length=0, trim_g=false, minimum_length=0, limit=0)

Parameters for sequence dereplication.

- `barcode_length > 0`: barcode at 5' end
- `barcode_length < 0`: barcode at 3' end
- `trim_g`: strip leading G nucleotides (RACE protocol artifact)
- `minimum_length`: discard reads shorter than this
- `limit`: process only first N reads (0 = all)
"""
struct DereplicateParams
    barcode_length::Int
    trim_g::Bool
    minimum_length::Int
    limit::Int
end

DereplicateParams(; barcode_length::Int=0, trim_g::Bool=false,
                    minimum_length::Int=0, limit::Int=0) =
    DereplicateParams(barcode_length, trim_g, minimum_length, limit)

"""
    DereplicatedRecord(name, sequence, barcode, group_size)

A dereplicated sequence with its barcode and the number of identical sequences found.
"""
struct DereplicatedRecord
    name::String
    sequence::String
    barcode::String
    group_size::Int
end

"""
    Dereplicator(params)

Callable struct that dereplicates FASTA/FASTQ sequences.

```julia
derep = Dereplicator(DereplicateParams(barcode_length=12, trim_g=true))
records, stats = derep("reads.fasta.gz")
```
"""
struct Dereplicator
    params::DereplicateParams
end

Dereplicator(; kwargs...) = Dereplicator(DereplicateParams(; kwargs...))

struct DereplicateStats
    total_reads::Int
    too_short::Int
    unique_sequences::Int
end

function (d::Dereplicator)(input_path::AbstractString)
    p = d.params
    sequences = Dict{String, Vector{FastaRecord}}()
    n_total = 0
    n_short = 0

    for rec in read_fasta(input_path)
        p.limit > 0 && n_total >= p.limit && break
        n_total += 1
        if length(rec.sequence) < p.minimum_length
            n_short += 1
            continue
        end
        push!(get!(sequences, rec.sequence, FastaRecord[]), rec)
    end

    results = DereplicatedRecord[]
    for (seq, recs) in sequences
        first_rec = first(recs)
        name = first(split(first_rec.name; limit=2))
        endswith(name, ';') && (name = name[1:end-1])

        barcode, unbarcoded = split_barcode(seq, p.barcode_length)
        p.trim_g && (unbarcoded = lstrip(unbarcoded, 'G'))

        push!(results, DereplicatedRecord(name, unbarcoded, barcode, length(recs)))
    end

    stats = DereplicateStats(n_total, n_short, length(results))
    (results, stats)
end

function split_barcode(sequence::AbstractString, barcode_length::Int)
    barcode_length == 0 && return ("", sequence)
    if barcode_length > 0
        bc = sequence[1:min(barcode_length, length(sequence))]
        seq = barcode_length < length(sequence) ? sequence[barcode_length+1:end] : ""
        (bc, seq)
    else
        cutpoint = max(1, length(sequence) + barcode_length + 1)
        bc = sequence[cutpoint:end]
        seq = sequence[1:cutpoint-1]
        (bc, seq)
    end
end

"""
    write_dereplicated(io, records, has_barcode)

Write dereplicated records in FASTA format with size annotation.
"""
function write_dereplicated(io::IO, records::Vector{DereplicatedRecord}, has_barcode::Bool)
    for rec in records
        if has_barcode
            println(io, ">$(rec.name);barcode=$(rec.barcode);size=$(rec.group_size);")
        else
            println(io, ">$(rec.name);size=$(rec.group_size);")
        end
        println(io, rec.sequence)
    end
end
