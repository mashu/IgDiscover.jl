# Gene renaming — assign canonical names to discovered sequences
# Maps each discovered sequence to its closest database gene name

"""
    rename_genes(discovered_fasta, database_fasta, output_fasta)

Rename discovered gene sequences based on closest match to original database.
Each discovered sequence gets the name of its closest database gene, with
allele suffixes (*01, *02, ...) assigned by cluster order.

This mirrors igdiscover's rename functionality for output normalization.
"""
function rename_genes(
    discovered_path::AbstractString,
    database_path::AbstractString,
    output_path::AbstractString,
)
    discovered = read_fasta(discovered_path)
    isempty(discovered) && (write_fasta(output_path, FastaRecord[]); return)

    database = read_fasta(database_path)
    db_seqs = Dict(r.name => r.sequence for r in database)

    # For each discovered sequence, find the closest database gene
    closest = Dict{String,Vector{Tuple{String,String,Int}}}()  # gene_base → [(name, seq, dist)]

    for rec in discovered
        best_name = ""
        best_dist = typemax(Int)

        for db_rec in database
            d = edit_distance(rec.sequence, db_rec.sequence; maxdiff=best_dist)
            if d < best_dist
                best_dist = d
                best_name = db_rec.name
            end
        end

        # Extract gene base (before '*')
        gene_base = occursin('*', best_name) ? first(split(best_name, '*')) : best_name
        push!(get!(closest, gene_base, Tuple{String,String,Int}[]),
              (rec.name, rec.sequence, best_dist))
    end

    # Assign allele numbers within each gene group, ordered by distance to database
    renamed = FastaRecord[]
    for (gene_base, alleles) in sort!(collect(closest); by=first)
        sort!(alleles; by=x -> x[3])  # sort by distance to database
        for (i, (orig_name, seq, dist)) in enumerate(alleles)
            allele_num = @sprintf("%02d", i)
            new_name = "$(gene_base)*$(allele_num)"
            push!(renamed, FastaRecord(new_name, seq))
        end
    end

    write_fasta(output_path, renamed)
    @info "Renamed $(length(renamed)) sequences to $output_path"
end
