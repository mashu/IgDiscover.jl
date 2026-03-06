# PCR bias correction: group by UMI (barcode) ± (pseudo-)CDR3, output one sequence per group
# Mirrors IgDiscover22 group step (https://gitlab.com/gkhlab/igdiscover22/-/blob/main/src/igdiscover/cli/group.py)

const MIN_CONSENSUS_SEQUENCES = 3
const GROUP_CONSENSUS_THRESHOLD = 0.501
const CDR3_LENGTH_TOLERANCE = 2

"""
    extract_barcode(record::FastaRecord, barcode_length::Int) -> (barcode::String, unbarcoded::FastaRecord)

Extract barcode from sequence. Positive length = 5' barcode (first N bases);
negative = 3' barcode (last N bases). Returns barcode string and record with remaining sequence.
"""
function extract_barcode(record::FastaRecord, barcode_length::Int)
    seq = record.sequence
    n = length(seq)
    abs_len = abs(barcode_length)
    abs_len <= n || return (seq, record)  # too short; treat whole as "unbarcoded"
    if barcode_length > 0
        barcode = seq[1:abs_len]
        unbarcoded_seq = seq[abs_len+1:end]
    else
        barcode = seq[n - abs_len + 1:n]
        unbarcoded_seq = seq[1:n - abs_len]
    end
    (barcode, FastaRecord(record.name, unbarcoded_seq))
end

"""
    trim_leading_g(seq::String) -> String

Remove leading run of G nucleotides (RACE protocol artifact).
"""
function trim_leading_g(seq::AbstractString)
    i = 1
    while i <= length(seq) && (seq[i] == 'G' || seq[i] == 'g')
        i += 1
    end
    i > 1 ? seq[i:end] : seq
end

"""
    pseudo_cdr3(sequence::AbstractString, start_from_end::Int, end_from_end::Int) -> String

Extract pseudo-CDR3 as slice from 3' end. start_from_end > end_from_end (e.g. 80 and 60).
Returns sequence[(len - start_from_end + 1):(len - end_from_end)].
"""
function pseudo_cdr3(sequence::AbstractString, start_from_end::Int, end_from_end::Int)
    L = length(sequence)
    start_idx = max(1, L - start_from_end + 1)
    end_idx = max(1, L - end_from_end)
    start_idx <= end_idx ? sequence[start_idx:end_idx] : ""
end

"""
    hamming_neighbors(s::String) -> Vector{String}

Return s and all strings at Hamming distance 1 (for CDR3 clustering).
"""
function hamming_neighbors(s::AbstractString)
    n = length(s)
    out = String[s]
    for i in 1:n
        for c in "ACGT"
            s[i] == c && continue
            push!(out, s[1:i-1] * c * s[i+1:end])
        end
    end
    out
end

"""
    cluster_by_cdr3(records_with_cdr3::Vector{Tuple{FastaRecord,String}}, length_tolerance::Int=2) -> Vector{Vector{Tuple{FastaRecord,String}}}

Single-linkage clustering: link (pseudo-)CDR3s with Hamming distance ≤ 1 and sequence length diff ≤ length_tolerance.
Splits CDR3 clusters by length where adjacent lengths differ by more than length_tolerance.
"""
function cluster_by_cdr3(records_with_cdr3::Vector{Tuple{FastaRecord,String}}, length_tolerance::Int = CDR3_LENGTH_TOLERANCE)
    isempty(records_with_cdr3) && return Vector{Vector{Tuple{FastaRecord,String}}}()
    length(records_with_cdr3) == 1 && return [records_with_cdr3]

    cdr3_to_records = Dict{String,Vector{Tuple{FastaRecord,String}}}()
    for (rec, cdr3) in records_with_cdr3
        push!(get!(cdr3_to_records, cdr3, Tuple{FastaRecord,String}[]), (rec, cdr3))
    end
    unique_cdr3s = sort(collect(keys(cdr3_to_records)))

    # Single-linkage: two CDR3s link if Hamming distance ≤ 1
    linked(s::String, t::String) = length(s) == length(t) && hamming_distance(s, t) <= 1
    cdr3_clusters = single_linkage(unique_cdr3s, linked)

    # Within each CDR3 cluster, split by sequence length (adjacent length diff > length_tolerance)
    out = Vector{Vector{Tuple{FastaRecord,String}}}()
    for cdr3_cluster in cdr3_clusters
        all_recs = vcat([cdr3_to_records[c] for c in cdr3_cluster]...)
        sort!(all_recs; by = x -> length(x[1].sequence))
        cluster = Tuple{FastaRecord,String}[]
        prev_len = nothing
        for (rec, cdr3) in all_recs
            L = length(rec.sequence)
            if prev_len !== nothing && L > prev_len + length_tolerance
                push!(out, cluster)
                cluster = Tuple{FastaRecord,String}[]
            end
            push!(cluster, (rec, cdr3))
            prev_len = L
        end
        !isempty(cluster) && push!(out, cluster)
    end
    out
end

"""
    simple_consensus_for_group(sequences::Vector{String}; threshold::Float64=0.501) -> String

Compute consensus from unaligned sequences (column-wise majority). Uses multialign + consensus_sequence.
Returns consensus; may contain 'N' if ambiguous.
"""
function simple_consensus_for_group(sequences::Vector{String}; threshold::Float64 = GROUP_CONSENSUS_THRESHOLD)
    isempty(sequences) && return ""
    length(sequences) == 1 && return sequences[1]
    seqs = Dict(string(i) => s for (i, s) in enumerate(sequences))
    aligned = multialign(seqs; program = "muscle-fast")
    consensus_sequence(collect(values(aligned)); threshold = threshold, ambiguous = 'N')
end

"""
    pick_group_representative(cluster::Vector{Tuple{FastaRecord,String}}, barcode::String,
                             use_consensus::Bool, consensus_counter::Ref{Int},
                             group_by_cdr3::Bool) -> FastaRecord

For a single group: 1 seq → that seq; 2 → random; 3+ → consensus (if use_consensus and no N) else random.
Output name has ;barcode=X;size=N (and ;cdr3=Y if group_by_cdr3).
"""
function pick_group_representative(cluster::Vector{Tuple{FastaRecord,String}}, barcode::String,
                                  use_consensus::Bool, consensus_counter::Ref{Int},
                                  group_by_cdr3::Bool)
    n = length(cluster)
    if n == 1
        rec, cdr3 = cluster[1]
        seq = rec.sequence
        name = first(split(rec.name, r"\s+"))
        cdr3_str = group_by_cdr3 ? cdr3 : ""
        return make_group_header_and_record(name, seq, barcode, n, cdr3_str)
    end
    if n < MIN_CONSENSUS_SEQUENCES || !use_consensus
        rec, cdr3 = cluster[rand(1:n)]
        name = first(split(rec.name, r"\s+"))
        cdr3_str = group_by_cdr3 ? cdr3 : ""
        return make_group_header_and_record(name, rec.sequence, barcode, n, cdr3_str)
    end
    # n >= 3 and use_consensus
    seqs = [r[1].sequence for r in cluster]
    cons = simple_consensus_for_group(seqs; threshold = GROUP_CONSENSUS_THRESHOLD)
    if occursin('N', cons)
        rec, cdr3 = cluster[1]
        name = first(split(rec.name, r"\s+"))
        cdr3_str = group_by_cdr3 ? cdr3 : ""
        return make_group_header_and_record(name, rec.sequence, barcode, n, cdr3_str)
    end
    consensus_counter[] += 1
    name = "consensus$(consensus_counter[])"
    cdr3_str = group_by_cdr3 ? first(cluster)[2] : ""  # use first CDR3 as label
    make_group_header_and_record(name, cons, barcode, n, cdr3_str)
end

function make_group_header_and_record(name::String, sequence::String, barcode::String, size::Int, cdr3::String)
    header = if isempty(cdr3)
        "$name;barcode=$barcode;size=$size;"
    else
        "$name;barcode=$barcode;cdr3=$cdr3;size=$size;"
    end
    FastaRecord(header, sequence)
end

"""
    group_reads(input_path::AbstractString, output_path::AbstractString, config::Config;
                limit::Int=0) -> String

Group sequences by barcode (± pseudo/real CDR3), output one sequence per group (PCR bias correction).
When config.barcode_length == 0, copies input to output and returns output_path.
Otherwise: extract barcode, optional trim-G, group by barcode and optionally by CDR3,
pick representative per group, write FASTA with headers ;barcode=...;size=... (and ;cdr3=... if grouping by CDR3).
Returns output_path.
"""
function group_reads(input_path::AbstractString, output_path::AbstractString, config::Config;
                     limit::Int = 0)
    if config.barcode_length == 0
        return input_path  # pipeline will use original reads
    end

    records = read_fasta(input_path; limit = limit)
    min_len = config.group_minimum_length
    bl = config.barcode_length
    trim_g = config.group_trim_g
    group_by_cdr3 = config.group_by_cdr3 in ("pseudo", "real")
    use_consensus = config.barcode_consensus

    # Collect by barcode: barcode -> list of (unbarcoded FastaRecord, optional cdr3)
    barcode_groups = Dict{String,Vector{FastaRecord}}()
    too_short = 0
    for rec in records
        length(rec.sequence) < min_len && (too_short += 1; continue)
        barcode, unbarcoded = extract_barcode(rec, bl)
        if trim_g
            unbarcoded = FastaRecord(unbarcoded.name, trim_leading_g(unbarcoded.sequence))
        end
        push!(get!(barcode_groups, barcode, FastaRecord[]), unbarcoded)
    end

    @info "Grouping: $(length(records)) reads, $too_short too short, $(length(barcode_groups)) unique barcodes"

    consensus_counter = Ref(0)
    out_records = FastaRecord[]

    for barcode in sort(collect(keys(barcode_groups)))
        seqs = barcode_groups[barcode]
        # Duplicate names check (as in Python)
        names = [first(split(r.name, r"\s+")) for r in seqs]
        length(names) != length(unique(names)) && error("Duplicate sequence names in input")

        if !group_by_cdr3
            clusters = [[(r, "") for r in seqs]]
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
            isempty(records_with_cdr3) && continue
            clusters = cluster_by_cdr3(records_with_cdr3)
        end

        for cluster in clusters
            rec_out = pick_group_representative(cluster, barcode, use_consensus, consensus_counter, group_by_cdr3)
            push!(out_records, rec_out)
        end
    end

    if endswith(output_path, ".gz")
        open(output_path, "w") do io
            stream = GzipCompressorStream(io)
            for r in out_records
                println(stream, ">", r.name)
                println(stream, r.sequence)
            end
            close(stream)
        end
    else
        write_fasta(output_path, out_records)
    end

    @info "Grouping wrote $(length(out_records)) sequences to $output_path"
    output_path
end
