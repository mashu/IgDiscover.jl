# Augment AIRR table with IgDiscover-specific columns
# Faithfully replicates Python igdiscover augment behavior

# ─── Safe column access helpers (avoid missing propagation) ───

col_str(df::DataFrame, col::Symbol, i::Int) =
    hasproperty(df, col) ? (let v = df[i, col]; ismissing(v) ? "" : string(v) end) : ""

col_int(df::DataFrame, col::Symbol, i::Int) =
    hasproperty(df, col) ? (let v = df[i, col]; ismissing(v) ? 0 : Int(v) end) : 0

col_float(df::DataFrame, col::Symbol, i::Int, default::Float64 = 0.0) =
    hasproperty(df, col) ? (let v = df[i, col]; ismissing(v) ? default : Float64(v) end) : default

"""
    parse_header(header::AbstractString) -> (name::String, count::Int, barcode::String)

Extract size= and barcode= fields from FASTA header. count=0 if not found.
"""
function parse_header(header::AbstractString)
    fields = split(header, ';')
    name = String(fields[1])
    count = 0
    barcode = ""
    for f in fields[2:end]
        isempty(f) && continue
        if occursin('=', f)
            k, v = split(f, '='; limit = 2)
            k == "size" && (count = parse(Int, v))
            k == "barcode" && (barcode = String(v))
        end
    end
    (name, count, barcode)
end

"""
    count_alignment_errors(germline_aln, sequence_aln) -> Int

Count mismatches between two aligned sequences (same length). Returns 0 for empty.
"""
function count_alignment_errors(germline_aln::AbstractString, sequence_aln::AbstractString)
    (isempty(germline_aln) || isempty(sequence_aln)) && return 0
    count(((g, s),) -> g != s, zip(germline_aln, sequence_aln))
end

"""
    alignment_coverage(germline_aln::AbstractString, db_len::Int) -> Float64

Coverage = (ungapped germline alignment length) / db_len * 100.
"""
function alignment_coverage(germline_aln::AbstractString, db_len::Int)
    (isempty(germline_aln) || db_len == 0) && return 0.0
    ungapped = count(!=('-'), germline_aln)
    100.0 * ungapped / db_len
end

"""
    query_position(gene_germline_start, gene_sequence_start,
                   germline_alignment, sequence_alignment, ref_pos) -> Int

Given a 0-based position on the reference, return the 0-based position on the query.
Returns -1 if the alignment doesn't cover that position.
"""
function query_position(germline_start::Int, sequence_start::Int,
                       germline_alignment::AbstractString,
                       sequence_alignment::AbstractString,
                       reference_position::Int)
    ref_pos = germline_start - 1   # 0-based
    query_pos = sequence_start - 1 # 0-based

    ref_pos == reference_position && return query_pos

    for (ref_c, query_c) in zip(germline_alignment, sequence_alignment)
        ref_c != '-' && (ref_pos += 1)
        query_c != '-' && (query_pos += 1)
        ref_pos == reference_position && return query_pos
    end
    -1  # not found
end

"""
    augment_table(airr_df::DataFrame, database_dir::AbstractString;
                 sequence_type::String="Ig") -> DataFrame

Add IgDiscover-specific columns to an AIRR-formatted IgBLAST table.
This must match Python igdiscover augment output exactly.
"""
function augment_table(airr_df::DataFrame, database_dir::AbstractString;
                      sequence_type::String = "Ig")
    db_v = read_fasta_dict(joinpath(database_dir, "V.fasta"))
    db_d = read_fasta_dict(joinpath(database_dir, "D.fasta"))
    db_j = read_fasta_dict(joinpath(database_dir, "J.fasta"))

    # Precompute CDR3 start positions for all V genes and CDR3 end for all J genes
    v_cdr3_starts = Dict{String, Dict{String, Int}}()  # locus → gene → position
    j_cdr3_ends = Dict{String, Dict{String, Int}}()
    for locus in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
        v_starts = Dict{String,Int}()
        for (name, seq) in db_v
            pos = cdr3_start_in_v(seq, locus)
            pos > 0 && (v_starts[name] = pos)
        end
        v_cdr3_starts[locus] = v_starts

        j_ends = Dict{String,Int}()
        for (name, seq) in db_j
            pos = cdr3_end_in_j(seq, locus)
            pos > 0 && (j_ends[name] = pos)
        end
        j_cdr3_ends[locus] = j_ends
    end

    n = nrow(airr_df)
    result = copy(airr_df)

    # Strip % prefix from gene calls (added by makeblastdb to avoid GenBank name mangling)
    for col in (:v_call, :d_call, :j_call)
        hasproperty(result, col) || continue
        result[!, col] = [lstrip(coalesce(v, ""), '%') for v in result[!, col]]
    end

    # Parse headers for count and barcode
    counts = zeros(Int, n)
    barcodes = fill("", n)
    for i in 1:n
        seqid = col_str(result, :sequence_id, i)
        name, cnt, bc = parse_header(seqid)
        result.sequence_id[i] = name
        counts[i] = cnt
        barcodes[i] = bc
    end
    result.count = counts
    result.barcode = barcodes

    # V_SHM and J_SHM from identity columns (Python: 100 - v_identity)
    result.V_SHM = [begin
        v = col_float(result, :v_identity, i, NaN)
        isnan(v) ? 0.0 : 100.0 - v
    end for i in 1:n]

    result.J_SHM = [begin
        v = col_float(result, :j_identity, i, NaN)
        isnan(v) ? 0.0 : 100.0 - v
    end for i in 1:n]

    # Errors
    result.V_errors = [begin
        g = col_str(result, :v_germline_alignment, i)
        s = col_str(result, :v_sequence_alignment, i)
        count_alignment_errors(g, s)
    end for i in 1:n]

    result.D_errors = [begin
        g = col_str(result, :d_germline_alignment, i)
        s = col_str(result, :d_sequence_alignment, i)
        count_alignment_errors(g, s)
    end for i in 1:n]

    result.J_errors = [begin
        g = col_str(result, :j_germline_alignment, i)
        s = col_str(result, :j_sequence_alignment, i)
        count_alignment_errors(g, s)
    end for i in 1:n]

    # Coverage
    result.V_covered = [begin
        g = col_str(result, :v_germline_alignment, i)
        vc = col_str(result, :v_call, i)
        dlen = haskey(db_v, vc) ? length(db_v[vc]) : 0
        alignment_coverage(g, dlen)
    end for i in 1:n]

    result.D_covered = [begin
        g = col_str(result, :d_germline_alignment, i)
        dc = col_str(result, :d_call, i)
        dlen = haskey(db_d, dc) ? length(db_d[dc]) : 0
        alignment_coverage(g, dlen)
    end for i in 1:n]

    result.J_covered = [begin
        g = col_str(result, :j_germline_alignment, i)
        jc = col_str(result, :j_call, i)
        dlen = haskey(db_j, jc) ? length(db_j[jc]) : 0
        alignment_coverage(g, dlen)
    end for i in 1:n]

    # V_CDR3_start — the CDR3 start offset within the V alignment (0-based, for discovery)
    result.V_CDR3_start = zeros(Int, n)
    # Also fill cdr3 / cdr3_aa columns from reference positions
    hasproperty(result, :cdr3) || (result.cdr3 = fill("", n))
    hasproperty(result, :cdr3_aa) || (result.cdr3_aa = fill("", n))

    for i in 1:n
        vc = col_str(result, :v_call, i)
        jc = col_str(result, :j_call, i)
        locus = col_str(result, :locus, i)
        seq = col_str(result, :sequence, i)

        (isempty(vc) || isempty(jc) || isempty(locus) || isempty(seq)) && continue

        # CDR3 start from V database
        v_starts_for_locus = get(v_cdr3_starts, locus, Dict{String,Int}())
        cdr3_ref_start = get(v_starts_for_locus, vc, 0)
        cdr3_ref_start == 0 && continue

        v_gs = col_int(result, :v_germline_start, i)
        v_ss = col_int(result, :v_sequence_start, i)
        v_ge = col_int(result, :v_germline_end, i)
        v_se = col_int(result, :v_sequence_end, i)
        v_galn = col_str(result, :v_germline_alignment, i)
        v_saln = col_str(result, :v_sequence_alignment, i)

        (v_gs == 0 || v_ss == 0) && continue

        cdr3_query_start = query_position(v_gs, v_ss, v_galn, v_saln, cdr3_ref_start)
        if cdr3_query_start < 0
            # Rescue: assume alignment continues without indels
            cdr3_query_start = v_se + (cdr3_ref_start - v_ge)
        end

        # V_CDR3_start = offset of CDR3 start within V sequence (0-based)
        result.V_CDR3_start[i] = cdr3_query_start - v_ss + 1

        # CDR3 end from J database
        j_ends_for_locus = get(j_cdr3_ends, locus, Dict{String,Int}())
        cdr3_ref_end = get(j_ends_for_locus, jc, 0)
        cdr3_ref_end == 0 && continue

        j_gs = col_int(result, :j_germline_start, i)
        j_ss = col_int(result, :j_sequence_start, i)
        j_galn = col_str(result, :j_germline_alignment, i)
        j_saln = col_str(result, :j_sequence_alignment, i)

        (j_gs == 0 || j_ss == 0) && continue

        cdr3_query_end = query_position(j_gs, j_ss, j_galn, j_saln, cdr3_ref_end)
        cdr3_query_end < 0 && continue

        # Extract CDR3 nucleotide sequence (0-based positions → 1-based Julia)
        if 1 <= cdr3_query_start + 1 <= cdr3_query_end <= length(seq)
            cdr3_nt = seq[cdr3_query_start+1 : cdr3_query_end]
            result.cdr3[i] = cdr3_nt
            result.cdr3_aa[i] = translate(cdr3_nt)
            result.cdr3_start[i] = cdr3_query_start + 1
            result.cdr3_end[i] = cdr3_query_end
        end
    end

    # Ensure V_nt column (ungapped V sequence alignment)
    if hasproperty(result, :v_sequence_alignment)
        result.V_nt = replace.(coalesce.(result.v_sequence_alignment, ""), "-" => "")
    else
        result.V_nt = fill("", n)
    end

    # Ensure d_support column
    hasproperty(result, :d_support) || (result.d_support = fill(Inf, n))

    # Ensure stop_codon column
    if !hasproperty(result, :stop_codon)
        result.stop_codon = fill("F", n)
    end

    result
end
