# Augment AIRR table with IgDiscover-specific columns

# ─── Safe column access helpers ───

col_str(df::DataFrame, col::Symbol, i::Int) =
    hasproperty(df, col) ? coalesce(df[i, col], "") : ""

col_int(df::DataFrame, col::Symbol, i::Int) =
    hasproperty(df, col) ? coalesce(Int(df[i, col]), 0) : 0

col_float(df::DataFrame, col::Symbol, i::Int, default::Float64=0.0) =
    hasproperty(df, col) ? Float64(coalesce(df[i, col], default)) : default

"""
    parse_header(header) -> (name, count, barcode)

Extract size= and barcode= fields from a FASTA header. count=0 if not found.
"""
function parse_header(header::AbstractString)
    fields = split(header, ';')
    name = String(fields[1])
    count = 0
    barcode = ""
    for f in fields[2:end]
        isempty(f) && continue
        if occursin('=', f)
            k, v = split(f, '='; limit=2)
            k == "size" && (count = parse(Int, v))
            k == "barcode" && (barcode = String(v))
        end
    end
    (name, count, barcode)
end

"""
    count_alignment_errors(germline_aln, sequence_aln) -> Int

Count mismatches between two aligned sequences (gaps count as mismatches).
"""
function count_alignment_errors(germline_aln::AbstractString, sequence_aln::AbstractString)
    (isempty(germline_aln) || isempty(sequence_aln)) && return 0
    count(((g, s),) -> g != s, zip(germline_aln, sequence_aln))
end

"""
    alignment_coverage(germline_aln, db_len) -> Float64

Coverage = (ungapped germline alignment length) / db_len * 100.
"""
function alignment_coverage(germline_aln::AbstractString, db_len::Int)
    (isempty(germline_aln) || db_len == 0) && return 0.0
    100.0 * count(!=('-'), germline_aln) / db_len
end

"""
    query_position(germline_start, sequence_start, germline_aln, sequence_aln, ref_pos) -> Int

Given a 0-based reference position, return the corresponding 0-based query position.
Returns -1 if the alignment does not cover that position.
"""
function query_position(
    germline_start::Int, sequence_start::Int,
    germline_aln::AbstractString, sequence_aln::AbstractString,
    reference_position::Int,
)
    ref_pos = germline_start - 1
    query_pos = sequence_start - 1
    ref_pos == reference_position && return query_pos
    for (ref_c, query_c) in zip(germline_aln, sequence_aln)
        ref_c != '-' && (ref_pos += 1)
        query_c != '-' && (query_pos += 1)
        ref_pos == reference_position && return query_pos
    end
    -1
end

# Build short_id → full_name mapping for IMGT-style headers (e.g. M99641|IGHV1-18*01|...).
# Registers: "part1_part2", "part2" (gene+allele), and gene-only so IgBLAST call styles all resolve.
function build_short_id_map(db::Dict{String,String})
    out = Dict{String,String}()
    for name in keys(db)
        parts = split(name, '|')
        get!(out, name, name)
        if length(parts) >= 2
            short = join(parts[1:2], "_")
            get!(out, short, name)
            gene_allele = parts[2]  # e.g. IGHV1-18*01, IGHJ1*01
            get!(out, gene_allele, name)
            if occursin('*', gene_allele)
                gene_only = first(split(gene_allele, '*'))
                get!(out, gene_only, name)  # e.g. IGHV1-18, IGHJ1
            end
        end
    end
    out
end

# From full IMGT header return short gene identifier (second pipe field: IGHV/IGHD/IGHJ style).
full_to_short(full_name::AbstractString) =
    (p = split(full_name, '|'); length(p) >= 2 ? p[2] : full_name)

"""
    augment_table(airr_df, database_dir; sequence_type) -> DataFrame

Add IgDiscover-specific columns to an AIRR-format IgBLAST table.
"""
function augment_table(
    airr_df::DataFrame,
    database_dir::AbstractString;
    sequence_type::String="Ig",
)
    db_v = read_fasta_dict(joinpath(database_dir, "V.fasta"))
    db_d = read_fasta_dict(joinpath(database_dir, "D.fasta"))
    db_j = read_fasta_dict(joinpath(database_dir, "J.fasta"))

    short_v = build_short_id_map(db_v)
    short_d = build_short_id_map(db_d)
    short_j = build_short_id_map(db_j)

    # Precompute CDR3 anchor positions
    loci = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
    v_cdr3_starts = Dict(locus => Dict{String,Int}() for locus in loci)
    j_cdr3_ends   = Dict(locus => Dict{String,Int}() for locus in loci)
    for locus in loci
        for (name, seq) in db_v
            pos = cdr3_start_in_v(seq, locus)
            pos > 0 && (v_cdr3_starts[locus][name] = pos)
        end
        for (name, seq) in db_j
            pos = cdr3_end_in_j(seq, locus)
            pos > 0 && (j_cdr3_ends[locus][name] = pos)
        end
    end

    n = nrow(airr_df)
    result = copy(airr_df)

    # Strip % prefix from gene calls
    for col in (:v_call, :d_call, :j_call)
        hasproperty(result, col) || continue
        result[!, col] = [lstrip(coalesce(v, ""), '%') for v in result[!, col]]
    end

    # Parse headers for count and barcode
    counts   = zeros(Int, n)
    barcodes = fill("", n)
    for i in 1:n
        name, cnt, bc = parse_header(col_str(result, :sequence_id, i))
        result.sequence_id[i] = name
        counts[i]   = cnt
        barcodes[i] = bc
    end
    result.count   = counts
    result.barcode = barcodes

    # SHM = 100 - identity
    result.V_SHM = [let v = col_float(result, :v_identity, i, NaN)
        isnan(v) ? 0.0 : 100.0 - v end for i in 1:n]
    result.J_SHM = [let v = col_float(result, :j_identity, i, NaN)
        isnan(v) ? 0.0 : 100.0 - v end for i in 1:n]

    # Alignment error counts
    result.V_errors = [count_alignment_errors(
        col_str(result, :v_germline_alignment, i), col_str(result, :v_sequence_alignment, i)) for i in 1:n]
    result.D_errors = [count_alignment_errors(
        col_str(result, :d_germline_alignment, i), col_str(result, :d_sequence_alignment, i)) for i in 1:n]
    result.J_errors = [count_alignment_errors(
        col_str(result, :j_germline_alignment, i), col_str(result, :j_sequence_alignment, i)) for i in 1:n]

    # Gene coverage
    resolve_v(vc) = get(short_v, vc, vc)
    resolve_d(dc) = get(short_d, dc, dc)
    resolve_j(jc) = get(short_j, jc, jc)

    result.V_covered = [alignment_coverage(
        col_str(result, :v_germline_alignment, i),
        get(db_v, resolve_v(col_str(result, :v_call, i)), "") |> length) for i in 1:n]
    result.D_covered = [alignment_coverage(
        col_str(result, :d_germline_alignment, i),
        get(db_d, resolve_d(col_str(result, :d_call, i)), "") |> length) for i in 1:n]
    result.J_covered = [alignment_coverage(
        col_str(result, :j_germline_alignment, i),
        get(db_j, resolve_j(col_str(result, :j_call, i)), "") |> length) for i in 1:n]

    # Ensure CDR3/stop columns exist as mutable String vectors
    for col in (:cdr3, :cdr3_aa)
        hasproperty(result, col) || (result[!, col] = fill("", n))
        result[!, col] = Vector{String}(coalesce.(result[!, col], ""))
    end
    for col in (:cdr3_start, :cdr3_end)
        hasproperty(result, col) || (result[!, col] = zeros(Int, n))
        result[!, col] = Vector{Int}(coalesce.(result[!, col], 0))
    end
    result.V_CDR3_start = zeros(Int, n)

    # CDR3 extraction from V/J anchor positions
    for i in 1:n
        vc     = col_str(result, :v_call, i)
        jc     = col_str(result, :j_call, i)
        locus  = col_str(result, :locus, i)
        seq    = col_str(result, :sequence, i)
        (isempty(vc) || isempty(jc) || isempty(locus) || isempty(seq)) && continue

        vc_full = resolve_v(vc)
        jc_full = resolve_j(jc)

        cdr3_ref_start = get(get(v_cdr3_starts, locus, Dict{String,Int}()), vc_full, 0)
        cdr3_ref_start == 0 && continue

        v_gs   = col_int(result, :v_germline_start, i)
        v_ss   = col_int(result, :v_sequence_start, i)
        v_ge   = col_int(result, :v_germline_end, i)
        v_se   = col_int(result, :v_sequence_end, i)
        v_galn = col_str(result, :v_germline_alignment, i)
        v_saln = col_str(result, :v_sequence_alignment, i)
        (v_gs == 0 || v_ss == 0) && continue

        cdr3_query_start = query_position(v_gs, v_ss, v_galn, v_saln, cdr3_ref_start)
        if cdr3_query_start < 0
            cdr3_query_start = v_se + (cdr3_ref_start - v_ge)
        end
        result.V_CDR3_start[i] = cdr3_query_start - v_ss + 1

        cdr3_ref_end = get(get(j_cdr3_ends, locus, Dict{String,Int}()), jc_full, 0)
        cdr3_ref_end == 0 && continue

        j_gs   = col_int(result, :j_germline_start, i)
        j_ss   = col_int(result, :j_sequence_start, i)
        j_galn = col_str(result, :j_germline_alignment, i)
        j_saln = col_str(result, :j_sequence_alignment, i)
        (j_gs == 0 || j_ss == 0) && continue

        cdr3_query_end = query_position(j_gs, j_ss, j_galn, j_saln, cdr3_ref_end)
        cdr3_query_end < 0 && continue

        if 1 <= cdr3_query_start + 1 <= cdr3_query_end <= length(seq)
            cdr3_nt = seq[cdr3_query_start+1:cdr3_query_end]
            result.cdr3[i]       = cdr3_nt
            result.cdr3_aa[i]    = translate(cdr3_nt)
            result.cdr3_start[i] = cdr3_query_start + 1
            result.cdr3_end[i]   = cdr3_query_end
        end
    end

    # V_nt: ungapped V sequence alignment
    result.V_nt = hasproperty(result, :v_sequence_alignment) ?
        replace.(coalesce.(result.v_sequence_alignment, ""), "-" => "") :
        fill("", n)

    # d_support: D-gene e-value
    hasproperty(result, :d_support) || (result.d_support = fill(Inf, n))

    # Normalize v_call, j_call, d_call to short identifiers (e.g. IGHV1-18*01, IGHD1-1*01, IGHJ1*01)
    for col in (:v_call, :d_call, :j_call)
        hasproperty(result, col) || continue
        resolve = col === :v_call ? resolve_v : (col === :d_call ? resolve_d : resolve_j)
        result[!, col] = [full_to_short(resolve(coalesce(result[i, col], ""))) for i in 1:n]
    end

    # Normalize stop_codon to "T"/"F"
    if !hasproperty(result, :stop_codon)
        result.stop_codon = fill("F", n)
    else
        result.stop_codon = [let s = lowercase(string(coalesce(x, "")))
            s == "true" ? "T" : (s == "false" ? "F" : string(x))
        end for x in result.stop_codon]
    end

    result
end
