# PCR bias correction: group by UMI (barcode) ± (pseudo-)CDR3, output one sequence per group

const MIN_CONSENSUS_SEQUENCES = 3
const GROUP_CONSENSUS_THRESHOLD = 0.501
const CDR3_LENGTH_TOLERANCE = 2

# ─── ConsensusCounter functor (replaces Ref{Int}) ───

"""
    ConsensusCounter()

Callable that generates unique consensus sequence names.
Each call increments the internal counter and returns "consensusN".
"""
mutable struct ConsensusCounter
    count::Int
    ConsensusCounter() = new(0)
end

function (cc::ConsensusCounter)()
    cc.count += 1
    "consensus$(cc.count)"
end

# ─── Helpers ───

"""
    extract_barcode(record, barcode_length) -> (barcode, unbarcoded_record)

Extract barcode from sequence. Positive length = 5' barcode; negative = 3' barcode.
"""
function extract_barcode(record::FastaRecord, barcode_length::Int)
    seq = record.sequence
    n = length(seq)
    abs_len = abs(barcode_length)
    abs_len <= n || return (seq, record)
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
"""
function pseudo_cdr3(sequence::AbstractString, start_from_end::Int, end_from_end::Int)
    L = length(sequence)
    start_idx = max(1, L - start_from_end + 1)
    end_idx = max(1, L - end_from_end)
    start_idx <= end_idx ? sequence[start_idx:end_idx] : ""
end

"""
    cluster_by_cdr3(records_with_cdr3, length_tolerance) -> Vector of clusters

Single-linkage clustering: link CDR3s with Hamming distance <= 1, then split by read length.
"""
function cluster_by_cdr3(
    records_with_cdr3::Vector{Tuple{FastaRecord,String}},
    length_tolerance::Int=CDR3_LENGTH_TOLERANCE,
)
    isempty(records_with_cdr3) && return Vector{Vector{Tuple{FastaRecord,String}}}()
    length(records_with_cdr3) == 1 && return [records_with_cdr3]

    cdr3_to_records = Dict{String,Vector{Tuple{FastaRecord,String}}}()
    for (rec, cdr3) in records_with_cdr3
        push!(get!(cdr3_to_records, cdr3, Tuple{FastaRecord,String}[]), (rec, cdr3))
    end
    unique_cdr3s = sort!(collect(keys(cdr3_to_records)))

    linked(s::String, t::String) = length(s) == length(t) && hamming_distance(s, t) <= 1
    cdr3_clusters = single_linkage(unique_cdr3s, linked)

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

function make_group_record(name::String, sequence::String, barcode::String, size::Int, cdr3::String)
    header = isempty(cdr3) ?
        "$name;barcode=$barcode;size=$size;" :
        "$name;barcode=$barcode;cdr3=$cdr3;size=$size;"
    FastaRecord(header, sequence)
end

"""
    pick_group_representative(cluster, barcode, use_consensus, counter, group_by_cdr3)

Choose a representative: 1-2 sequences -> pick first; 3+ -> consensus if unambiguous.
"""
function pick_group_representative(
    cluster::Vector{Tuple{FastaRecord,String}},
    barcode::String,
    use_consensus::Bool,
    counter::ConsensusCounter,
    group_by_cdr3::Bool,
    program::String="muscle-fast",
)
    n = length(cluster)
    cdr3_label = group_by_cdr3 ? cluster[1][2] : ""

    if n < MIN_CONSENSUS_SEQUENCES || !use_consensus
        rec = cluster[1][1]
        return make_group_record(first(split(rec.name, r"\s+")), rec.sequence, barcode, n, cdr3_label)
    end

    seqs = [r[1].sequence for r in cluster]
    cons = simple_consensus_for_group(seqs; program=program, threshold=GROUP_CONSENSUS_THRESHOLD)

    if occursin('N', cons)
        rec = cluster[1][1]
        return make_group_record(first(split(rec.name, r"\s+")), rec.sequence, barcode, n, cdr3_label)
    end

    make_group_record(counter(), cons, barcode, n, cdr3_label)
end

function simple_consensus_for_group(sequences::Vector{String}; program::String="muscle-fast", threshold::Float64=GROUP_CONSENSUS_THRESHOLD)
    isempty(sequences) && return ""
    length(sequences) == 1 && return sequences[1]
    seqs = Dict(string(i) => s for (i, s) in enumerate(sequences))
    aligned = multialign(seqs; program=program)
    consensus_sequence(collect(values(aligned)); threshold=threshold, ambiguous='N')
end

"""
    group_reads(input_path, output_path, config; limit) -> String

Group sequences by barcode (± pseudo/real CDR3) for PCR bias correction.
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

    counter = ConsensusCounter()
    out_records = FastaRecord[]

    for barcode in sort!(collect(keys(barcode_groups)))
        seqs = barcode_groups[barcode]

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
                cluster, barcode, config.barcode_consensus, counter, group_by_cdr3, config.multialign_program))
        end
    end

    write_fasta_gz(output_path, out_records)
    @info "Grouping wrote $(length(out_records)) sequences to $output_path"
    output_path
end

"""
    trim_reads_race_g(input_path, output_path; limit) -> String

Trim leading G nucleotides from reads (RACE protocol artifact).
"""
function trim_reads_race_g(
    input_path::AbstractString,
    output_path::AbstractString;
    limit::Int=0,
)
    isfile(output_path) && return output_path

    records = read_fasta(input_path; limit=limit)
    trimmed = [FastaRecord(r.name, trim_leading_g(r.sequence)) for r in records]
    filter!(r -> !isempty(r.sequence), trimmed)

    write_fasta_gz(output_path, trimmed)
    @info "Trimmed leading G from $(length(records)) reads -> $(length(trimmed)) in $output_path"
    output_path
end
