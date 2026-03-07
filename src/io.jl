# I/O utilities for FASTA files and assignment tables

"""
    FastaRecord

A single FASTA record with `name` (header/identifier) and `sequence` (nucleotide or amino acid string).
"""
struct FastaRecord
    name::String
    sequence::String
end

# Minimal IO wrapper that prepends one byte so we can peek and still pass stream to reader
mutable struct PrependIO <: IO
    first::Union{UInt8,Nothing}
    rest::IO
end
Base.read(io::PrependIO, ::Type{UInt8}) = io.first !== nothing ? (b = io.first; io.first = nothing; b) : read(io.rest, UInt8)
Base.eof(io::PrependIO) = io.first === nothing && eof(io.rest)
Base.close(io::PrependIO) = close(io.rest)
Base.isopen(io::PrependIO) = isopen(io.rest)
Base.bytesavailable(io::PrependIO) = (io.first !== nothing ? 1 : 0) + bytesavailable(io.rest)

"""
    read_fasta(path; limit=0) -> Vector{FastaRecord}

Read records from a FASTA or FASTQ file (plain or gzipped). Format is auto-detected
by the first character ('>' = FASTA, '@' = FASTQ). Sequences are uppercased.
If `limit > 0`, stop after that many records.
"""
function read_fasta(path::AbstractString; limit::Int=0)
    records = FastaRecord[]
    raw = endswith(path, ".gz") ? GzipDecompressorStream(open(path)) : open(path)
    first_chunk = read(raw, 1)
    isempty(first_chunk) && (close(raw); return records)
    first_byte = only(first_chunk)
    stream = PrependIO(first_byte, raw)
    if first_byte == UInt8('>')
        reader = FASTA.Reader(stream)
        for record in reader
            name = FASTA.identifier(record)
            seq = uppercase(String(FASTA.sequence(record)))
            push!(records, FastaRecord(name, seq))
            limit > 0 && length(records) >= limit && break
        end
    elseif first_byte == UInt8('@')
        reader = FASTQ.Reader(stream)
        for record in reader
            name = first(split(FASTQ.identifier(record), r"\s+"))
            seq = uppercase(String(FASTQ.sequence(record)))
            push!(records, FastaRecord(name, seq))
            limit > 0 && length(records) >= limit && break
        end
    else
        close(stream)
        error("Not FASTA (>) or FASTQ (@): first byte '$(Char(first_byte))' (0x$(string(first_byte, base=16))) in $path")
    end
    close(reader)
    records
end

"""
    read_fasta_dict(path) -> Dict{String,String}

Read FASTA into name → sequence dictionary.
"""
function read_fasta_dict(path::AbstractString)
    Dict(r.name => r.sequence for r in read_fasta(path))
end

# ─── IMGT database sanitization ───

"""
    allele_name_from_header(name) -> String

Extract allele identifier from IMGT-style header.
Examples:
  "J00256|IGHJ1*01|Homo sapiens|F|..."  → "IGHJ1*01"
  "IGHV1-18*01"                          → "IGHV1-18*01"
"""
allele_name_from_header(name::AbstractString) =
    (p = split(name, '|'); length(p) >= 2 ? String(p[2]) : String(name))

"""
    sanitize_imgt_sequence(seq) -> String

Remove IMGT gap characters (dots) from a nucleotide sequence and uppercase.
IMGT databases use '.' to represent gaps in their numbering system;
these must be removed before use with IgBLAST or alignment tools.
"""
sanitize_imgt_sequence(seq::AbstractString) = uppercase(replace(seq, "." => ""))

"""
    sanitize_imgt_record(record::FastaRecord) -> FastaRecord

Clean an IMGT FASTA record: extract allele name from header and remove dots from sequence.
"""
sanitize_imgt_record(record::FastaRecord) =
    FastaRecord(allele_name_from_header(record.name), sanitize_imgt_sequence(record.sequence))

"""
    write_sanitized_imgt(path_in, path_out)

Read IMGT FASTA from path_in and write to path_out with:
  - Headers replaced by allele names (e.g. IGHV1-18*01, IGHD3-3*01, IGHJ1*01)
  - IMGT gap dots removed from sequences
This ensures IgBLAST and downstream tools see clean identifiers and ungapped sequences.
"""
function write_sanitized_imgt(path_in::AbstractString, path_out::AbstractString)
    records = read_fasta(path_in)
    cleaned = sanitize_imgt_record.(records)
    # Remove empty sequences (can occur if IMGT entry is all gaps)
    filter!(r -> !isempty(r.sequence), cleaned)
    write_fasta(path_out, cleaned)
end

# ─── FASTA writing ───

"""
    write_fasta(path, records)

Write FASTA records to a plain-text file.
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
    write_fasta(path, [FastaRecord(k, v) for (k, v) in sort!(collect(d); by=first)])
end

"""
    write_fasta_gz(path, records)

Write FASTA records to a gzip-compressed file. Uses plain write_fasta for non-.gz paths.
"""
function write_fasta_gz(path::AbstractString, records::Vector{FastaRecord})
    if endswith(path, ".gz")
        open(path, "w") do io
            stream = GzipCompressorStream(io)
            for r in records
                println(stream, ">", r.name)
                println(stream, r.sequence)
            end
            close(stream)
        end
    else
        write_fasta(path, records)
    end
end

# ─── Assignment tables (TSV) ───

"""
    read_assignments(path) -> DataFrame

Read an assignment table (TSV, optionally gzipped).
"""
function read_assignments(path::AbstractString)
    df = CSV.read(path, DataFrame;
        delim='\t', missingstring=["", "NA"], ntasks=1, pool=false)
    for col in (:v_call, :d_call, :j_call, :locus, :cdr3, :sequence_id,
                :barcode, :sequence, :stop_codon)
        hasproperty(df, col) || continue
        df[!, col] = coalesce.(df[!, col], "")
    end
    df
end

"""
    write_table(path, df)
"""
function write_table(path::AbstractString, df::DataFrame)
    CSV.write(path, df; delim='\t')
end

"""
    write_table_gz(path, df)

Write DataFrame as gzip-compressed TSV.
"""
function write_table_gz(path::AbstractString, df::DataFrame)
    CSV.write(path, df; delim='\t', compress=true)
end

"""
    validate_fasta(path)

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
