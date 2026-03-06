# IgBLAST wrapper — calls external igblastn and makeblastdb
# Uses same external programs as Python igdiscover.
# Consider IgBLAST.jl (https://github.com/mashu/IgBLAST.jl) for binary management.

"""
    make_blastdb(fasta_path, db_name; prefix="%")

Create a BLAST database. The `%` prefix avoids IgBLAST GenBank-style name mangling.
"""
function make_blastdb(fasta_path::AbstractString, db_name::AbstractString;
                     prefix::String = "%")
    records = read_fasta(fasta_path)
    isempty(records) && error("FASTA file $fasta_path is empty")

    # Sanitize sequences for makeblastdb: only A,C,G,T,U,N allowed; IMGT gaps '.' etc. -> N
    function sanitize_nucl(s::AbstractString)
        buf = IOBuffer(sizehint = length(s))
        for c in uppercase(s)
            write(buf, c in "ACGTUN" ? c : 'N')
        end
        String(take!(buf))
    end
    db_fasta = db_name * ".fasta"
    open(db_fasta, "w") do io
        for r in records
            parts = split(r.name, '|')
            # Unique ID: accession_allele (IMGT has duplicate accessions across alleles)
            seq_id = length(parts) >= 2 ? join(parts[1:2], "_") : first(parts)
            println(io, ">", prefix, seq_id)
            println(io, sanitize_nucl(r.sequence))
        end
    end
    run(pipeline(`makeblastdb -parse_seqids -dbtype nucl -in $db_fasta -out $db_name`;
                 stdout=devnull, stderr=devnull))
end

"""
    make_vdj_blastdb(blastdb_dir, database_dir)

Build V, D, J BLAST databases in `blastdb_dir` from FASTA files in `database_dir`.
"""
function make_vdj_blastdb(blastdb_dir::AbstractString, database_dir::AbstractString)
    for gene in ("V", "D", "J")
        make_blastdb(joinpath(database_dir, gene * ".fasta"),
                     joinpath(blastdb_dir, gene))
    end
end

# ─── Two dispatch variants for penalty ───

"""
    run_igblast_chunk(sequences, blastdb_dir; species, sequence_type) -> String

Run igblastn on sequences (no mismatch penalty override).
"""
function run_igblast_chunk(sequences::Vector{FastaRecord},
                          blastdb_dir::AbstractString;
                          species::String = "",
                          sequence_type::String = "Ig")
    run_igblast_impl(sequences, blastdb_dir, species, sequence_type, String[])
end

"""
    run_igblast_chunk(sequences, blastdb_dir, penalty; species, sequence_type) -> String

Run igblastn on sequences with a mismatch penalty.
"""
function run_igblast_chunk(sequences::Vector{FastaRecord},
                          blastdb_dir::AbstractString,
                          penalty::Int;
                          species::String = "",
                          sequence_type::String = "Ig")
    run_igblast_impl(sequences, blastdb_dir, species, sequence_type, ["-penalty", string(penalty)])
end

function run_igblast_impl(sequences::Vector{FastaRecord},
                         blastdb_dir::String, species::String,
                         sequence_type::String, extra_args::Vector{String})
    sequence_type in ("Ig", "TCR") || error("sequence_type must be \"Ig\" or \"TCR\"")

    # Build argument list
    args = String[]
    for gene in ("V", "D", "J")
        append!(args, ["-germline_db_$(gene)", joinpath(blastdb_dir, gene)])
    end

    # Empty aux file suppresses warnings
    aux_path = tempname() * ".aux"
    touch(aux_path)
    append!(args, ["-auxiliary_data", aux_path])
    append!(args, extra_args)

    if !isempty(species)
        append!(args, ["-organism", species])
    end
    append!(args, ["-ig_seqtype", sequence_type,
                   "-num_threads", "1",
                   "-domain_system", "imgt",
                   "-num_alignments_V", "1",
                   "-num_alignments_D", "1",
                   "-num_alignments_J", "1",
                   "-outfmt", "19",
                   "-query", "-"])

    fasta_input = join(">$(r.name)\n$(r.sequence)\n" for r in sequences)

    outfile = tempname()
    args_with_out = vcat(args, ["-out", outfile])

    run(pipeline(Cmd(vcat("igblastn", args_with_out));
                 stdin=IOBuffer(fasta_input), stdout=devnull, stderr=devnull))
    result = read(outfile, String)
    rm(outfile; force = true)
    rm(aux_path; force = true)
    result
end

"""
    run_igblast_on_fasta(database_dir, input_fasta, output_tsv_gz;
                        sequence_type, species, penalty, threads, limit, chunksize)

Run IgBLAST on a FASTA file, write AIRR TSV (gzipped) output.
Main entry point for the IgBLAST step of the pipeline.
When `limit > 0`, only the first `limit` reads are processed (for testing).
"""
function run_igblast_on_fasta(database_dir::AbstractString,
                             input_fasta::AbstractString,
                             output_path::AbstractString;
                             sequence_type::String = "Ig",
                             species::String = "",
                             penalty::Int = 0,
                             threads::Int = Sys.CPU_THREADS,
                             limit::Int = 0,
                             chunksize::Int = 1000)
    sequences = read_fasta(input_fasta; limit = limit)
    limit_str = limit > 0 ? " (limit=$limit)" : ""
    @info "Running IgBLAST on $(length(sequences)) sequences with $threads threads$limit_str"

    blastdb_dir = mktempdir()
    make_vdj_blastdb(blastdb_dir, database_dir)

    chunks = [sequences[i:min(i + chunksize - 1, end)]
              for i in 1:chunksize:length(sequences)]

    results = Vector{String}(undef, length(chunks))

    # Parallel execution via Threads
    Threads.@threads for idx in eachindex(chunks)
        results[idx] = if penalty > 0
            run_igblast_chunk(chunks[idx], blastdb_dir, penalty;
                             species = species, sequence_type = sequence_type)
        else
            run_igblast_chunk(chunks[idx], blastdb_dir;
                             species = species, sequence_type = sequence_type)
        end
    end

    # Combine: keep header from first, skip header in rest
    io = open(output_path, "w")
    stream = GzipCompressorStream(io)
    for (i, result) in enumerate(results)
        for (j, line) in enumerate(eachline(IOBuffer(result)))
            (i > 1 && j == 1) && continue  # skip header of subsequent chunks
            isempty(line) && continue
            println(stream, line)
        end
    end
    close(stream)
    close(io)
    rm(blastdb_dir; recursive = true, force = true)

    @info "IgBLAST complete, output written to $output_path"
end
