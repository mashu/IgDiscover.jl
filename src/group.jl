# PCR bias correction: group by UMI (barcode) ± (pseudo-)CDR3, output one sequence per group
# Mirrors IgDiscover22 group step

const MIN_CONSENSUS_SEQUENCES = 3
const GROUP_CONSENSUS_THRESHOLD = 0.501
const CDR3_LENGTH_TOLERANCE = 2

"""
    extract_barcode(record, barcode_length) -> (barcode, unbarcoded_record)

Extract barcode from sequence. Positive length = 5' barcode; negative = 3' barcode.
Returns the barcode string and a record with the remaining sequence.
"""
function extract_barcode(record::FastaRecord, barcode_length::Int)
    seq = record.sequence
    n = length(seq)
    abs_len = abs(barcode_length)
    abs_len <= n || return (seq, record)  # too short; treat whole sequence as barcode
    if barcode_length > 0
        (seq[1:abs_len], FastaRecord(record.name, seq[abs_len+1:end]))
    else
        (seq[n-abs_len+1:n], FastaRecord(record.name, seq[1:n-abs_len]))
    end
end

"""
    trim_leading_g(seq) -> String

Remove leading run of G nucleotides (RACE protocol artifact).
"""
function trim_leading_g(seq::AbstractString)
    i = findfirst(c -> c != 'G' && c != 'g', seq)
    i === nothing ? "" : (i > 1 ? seq[i:end] : seq)
end

"""
    pseudo_cdr3(sequence, start_from_end, end_from_end) -> String

Extract pseudo-CDR3 as a slice from the 3' end.
Returns sequence[(len - start_from_end + 1):(len - end_from_end)].
"""
function pseudo_cdr3(sequence::AbstractString, start_from_end::Int, end_from_end::Int)
    L = length(sequence)
    start_idx = max(1, L - start_from_end + 1)
    end_idx = max(1, L - end_from_end)
    start_idx <= end_idx ? sequence[start_idx:end_idx] : ""
end

"""
    cluster_by_cdr3(records_with_cdr3, length_tolerance) -> Vector of clusters

Single-linkage clustering: link (pseudo-)CDR3s with Hamming distance ≤ 1 and
sequence length difference ≤ length_tolerance.
Each cluster is a Vector of (FastaRecord, cdr3_string) tuples.
"""
function cluster_by_cdr3(
    records_with_cdr3::Vector{Tuple{FastaRecord,String}},
    length_tolerance::Int=CDR3_LENGTH_TOLERANCE,
)
    isempty(records_with_cdr3) && return Vector{Vector{Tuple{FastaRecord,String}}}()
    length(records_with_cdr3) == 1 && return [records_with_cdr3]

    # Group records by CDR3 string
    cdr3_to_records = Dict{String,Vector{Tuple{FastaRecord,String}}}()
    for (rec, cdr3) in records_with_cdr3
        push!(get!(cdr3_to_records, cdr3, Tuple{FastaRecord,String}[]), (rec, cdr3))
    end
    unique_cdr3s = sort!(collect(keys(cdr3_to_records)))

    # Single-linkage on CDR3s: link if same length and Hamming distance ≤ 1
    linked(s::String, t::String) = length(s) == length(t) && hamming_distance(s, t) <= 1
    cdr3_clusters = single_linkage(unique_cdr3s, linked)

    # Within each CDR3 cluster, further split by read length (length jump > tolerance)
    out = Vector{Vector{Tuple{FastaRecord,String}}}()
    for cdr3_cluster in cdr3_clusters
        all_recs = reduce(vcat, cdr3_to_records[c] for c in cdr3_cluster)
        sort!(all_recs; by=x -> length(x[1].sequence))
        current_cluster = Tuple{FastaRecord,String}[]
        prev_len = -1
        for item in all_recs
            L = length(item[1].sequence)
            if prev_len >= 0 && L > prev_len + length_tolerance
                push!(out, current_cluster)
                current_cluster = Tuple{FastaRecord,String}[]
            end
            push!(current_cluster, item)
            prev_len = L
        end
        isempty(current_cluster) || push!(out, current_cluster)
    end
    out
end

"""
    make_group_header(name, sequence, barcode, size, cdr3) -> FastaRecord

Construct a FASTA record with IgDiscover group metadata in the header.
"""
function make_group_record(name::String, sequence::String, barcode::String, size::Int, cdr3::String)
    header = isempty(cdr3) ?
        "$name;barcode=$barcode;size=$size;" :
        "$name;barcode=$barcode;cdr3=$cdr3;size=$size;"
    FastaRecord(header, sequence)
end

"""
    pick_group_representative(cluster, barcode, use_consensus, consensus_counter, group_by_cdr3)
        -> FastaRecord

Choose a representative sequence for a barcode group:
- 1 sequence  → use it directly
- 2 sequences → pick first
- 3+          → compute consensus (if use_consensus and no ambiguous bases), else pick first
"""
function pick_group_representative(
    cluster::Vector{Tuple{FastaRecord,String}},
    barcode::String,
    use_consensus::Bool,
    consensus_counter::Ref{Int},
    group_by_cdr3::Bool,
)
    n = length(cluster)
    cdr3_label = group_by_cdr3 ? cluster[1][2] : ""

    if n < MIN_CONSENSUS_SEQUENCES || !use_consensus
        rec = cluster[1][1]
        return make_group_record(first(split(rec.name, r"\s+")), rec.sequence, barcode, n, cdr3_label)
    end

    seqs = [r[1].sequence for r in cluster]
    cons = simple_consensus_for_group(seqs; threshold=GROUP_CONSENSUS_THRESHOLD)

    if occursin('N', cons)
        # Consensus is ambiguous; fall back to first sequence
        rec = cluster[1][1]
        return make_group_record(first(split(rec.name, r"\s+")), rec.sequence, barcode, n, cdr3_label)
    end

    consensus_counter[] += 1
    make_group_record("consensus$(consensus_counter[])", cons, barcode, n, cdr3_label)
end

"""
    simple_consensus_for_group(sequences; threshold) -> String

Column-wise majority consensus from unaligned sequences via MUSCLE alignment.
May contain 'N' at ambiguous positions.
"""
function simple_consensus_for_group(sequences::Vector{String}; threshold::Float64=GROUP_CONSENSUS_THRESHOLD)
    isempty(sequences) && return ""
    length(sequences) == 1 && return sequences[1]
    seqs = Dict(string(i) => s for (i, s) in enumerate(sequences))
    aligned = multialign(seqs; program="muscle-fast")
    consensus_sequence(collect(values(aligned)); threshold=threshold, ambiguous='N')
end

"""
    group_reads(input_path, output_path, config; limit) -> String

Group sequences by barcode (± pseudo/real CDR3) for PCR bias correction.
When config.barcode_length == 0, returns input_path unchanged.
Otherwise writes one representative per group to output_path and returns output_path.
"""
function group_reads(
    input_path::AbstractString,
    output_path::AbstractString,
    config::Config;
    limit::Int=0,
)
    config.barcode_length == 0 && return input_path

    records = read_fasta(input_path; limit=limit)
    group_by_cdr3 = config.group_by_cdr3 in ("pseudo", "real")

    # Partition reads by barcode
    barcode_groups = Dict{String,Vector{FastaRecord}}()
    too_short = 0
    for rec in records
        length(rec.sequence) < config.group_minimum_length && (too_short += 1; continue)
        barcode, unbarcoded = extract_barcode(rec, config.barcode_length)
        if config.group_trim_g
            unbarcoded = FastaRecord(unbarcoded.name, trim_leading_g(unbarcoded.sequence))
        end
        push!(get!(barcode_groups, barcode, FastaRecord[]), unbarcoded)
    end
    @info "Grouping: $(length(records)) reads, $too_short too short, $(length(barcode_groups)) unique barcodes"

    consensus_counter = Ref(0)
    out_records = FastaRecord[]

    for barcode in sort!(collect(keys(barcode_groups)))
        seqs = barcode_groups[barcode]
        names = [first(split(r.name, r"\s+")) for r in seqs]
        length(names) == length(unique(names)) || error("Duplicate sequence names in input")

        clusters = if !group_by_cdr3
            [[(r, "") for r in seqs]]
        else
            records_with_cdr3 = Tuple{FastaRecord,String}[]
            for r in seqs
                cdr3 = if config.group_by_cdr3 == "real"
                    st, en = find_cdr3(r.sequence, "IGH")
                    st > 0 ? r.sequence[st:en] : continue
                else
                    pseudo_cdr3(r.sequence, config.pseudo_cdr3_start, config.pseudo_cdr3_end)
                end
                push!(records_with_cdr3, (r, cdr3))
            end
            isempty(records_with_cdr3) ? continue : cluster_by_cdr3(records_with_cdr3)
        end

        for cluster in clusters
            push!(out_records, pick_group_representative(
                cluster, barcode, config.barcode_consensus, consensus_counter, group_by_cdr3))
        end
    end

    write_fasta_gz(output_path, out_records)
    @info "Grouping wrote $(length(out_records)) sequences to $output_path"
    output_path
end
