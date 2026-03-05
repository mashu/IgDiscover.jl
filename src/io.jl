# I/O utilities for FASTA files and assignment tables

struct FastaRecord
    name::String
    sequence::String
end

"""
    read_fasta(path::AbstractString) -> Vector{FastaRecord}

Read all records from a FASTA file (plain or gzipped). Sequences are uppercased.
"""
function read_fasta(path::AbstractString)
    records = FastaRecord[]
    reader = if endswith(path, ".gz")
        FASTA.Reader(GzipDecompressorStream(open(path)))
    else
        open(FASTA.Reader, path)
    end
    for record in reader
        name = first(split(FASTA.identifier(record), r"\s+"))
        seq = uppercase(String(FASTA.sequence(record)))
        push!(records, FastaRecord(name, seq))
    end
    close(reader)
    records
end

"""
    read_fasta_dict(path::AbstractString) -> Dict{String,String}

Read FASTA into name→sequence dictionary.
"""
function read_fasta_dict(path::AbstractString)
    Dict(r.name => r.sequence for r in read_fasta(path))
end

"""
    write_fasta(path::AbstractString, records::Vector{FastaRecord})

Write FASTA records to file.
"""
function write_fasta(path::AbstractString, records::Vector{FastaRecord})
    open(path, "w") do io
        for r in records
            println(io, ">", r.name)
            println(io, r.sequence)
        end
    end
end

function write_fasta(path::AbstractString, d::Dict{String,String})
    records = [FastaRecord(k, v) for (k, v) in sort(collect(d); by = first)]
    write_fasta(path, records)
end

"""
    read_assignments(path::AbstractString) -> DataFrame

Read assignment table (TSV, optionally gzipped). CSV reads .gz natively.
String columns get empty string for missing.
"""
function read_assignments(path::AbstractString)
    df = CSV.read(path, DataFrame; delim = '\t', missingstring = ["", "NA"],
                  ntasks = 1, pool = false)

    string_cols = [:v_call, :d_call, :j_call, :locus, :cdr3, :sequence_id,
                   :barcode, :sequence, :stop_codon]
    for col in string_cols
        hasproperty(df, col) || continue
        df[!, col] = coalesce.(df[!, col], "")
    end
    df
end

"""
    write_table(path::AbstractString, df::DataFrame)
"""
function write_table(path::AbstractString, df::DataFrame)
    CSV.write(path, df; delim = '\t')
end

"""
    write_table_gzipped(path::AbstractString, df::DataFrame)

Write DataFrame as TSV with gzip compression. CSV handles compression natively.
"""
function write_table_gzipped(path::AbstractString, df::DataFrame)
    CSV.write(path, df; delim = '\t', compress = true)
end

"""
    validate_fasta(path::AbstractString)

Validate: no empty records, no duplicate names, no duplicate sequences.
"""
function validate_fasta(path::AbstractString)
    records = read_fasta(path)
    names = Set{String}()
    sequences = Dict{String,String}()
    for r in records
        isempty(r.sequence) && error("Record '$(r.name)' is empty in $path")
        r.name in names && error("Duplicate name '$(r.name)' in $path")
        push!(names, r.name)
        haskey(sequences, r.sequence) &&
            error("Records '$(r.name)' and '$(sequences[r.sequence])' are identical in $path")
        sequences[r.sequence] = r.name
    end
end
