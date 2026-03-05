# CDR3 detection — ported from igdiscover/species.py
# Detects CDR3 start within V genes and CDR3 end within J genes

# ─── CDR3 start in V genes ───

const CDR3_START_VH = r"[FY][FHVWY]C([ADEGIKMNRSTV\*]|$)"
const CDR3_START_VH_ALT = r"C(.[RK])"
const CDR3_START_VK = r"[FSVY][CFHNVY][CDFGLSW]"
const CDR3_START_VL = r"[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]"
const CDR3_START_TR = r"[YFH]C"

"""
    cdr3_start_in_v(nt::AbstractString, locus::AbstractString) -> Int

Find CDR3 start position (0-based nucleotide offset) within a V gene sequence.
Returns 0 if not found.
"""
function cdr3_start_in_v(nt::AbstractString, locus::AbstractString)
    aa = translate(nt)
    isempty(aa) && return 0

    if locus == "IGH"
        return cdr3_start_heavy(aa)
    elseif locus in ("IGK", "IGL", "TRG", "TRD")
        regex = if locus == "IGK"
            CDR3_START_VK
        elseif locus == "IGL"
            CDR3_START_VL
        else
            CDR3_START_TR
        end
        # Search in the last 15 aa
        head_len = max(0, length(aa) - 15)
        tail = SubString(aa, head_len + 1)
        m = match(regex, tail)
        m === nothing && return 0
        return 3 * (head_len + m.offset + length(m.match) - 1)
    elseif locus in ("TRA", "TRB")
        head_len = max(0, length(aa) - 8)
        tail = SubString(aa, head_len + 1)
        pos = findfirst('C', tail)
        pos === nothing && return 0
        return 3 * (head_len + pos)
    end
    0
end

function cdr3_start_heavy(aa::AbstractString)
    head_len = max(0, length(aa) - 15)
    tail = SubString(aa, head_len + 1)

    # Primary regex: [FY][FHVWY]C followed by start of CDR3
    m = match(CDR3_START_VH, tail)
    if m === nothing
        m = match(CDR3_START_VH_ALT, tail)
    end
    m === nothing && return 0

    # Find the start of the captured group (the CDR3 start residue)
    # In the primary regex, CDR3 starts after the 'C' (3 chars from match start)
    # We need the position of the capture group
    if length(m.captures) >= 1 && m.captures[1] !== nothing
        # Capture group offset within the match
        cap_offset = m.offsets[1]
        return 3 * (head_len + cap_offset - 1)
    end
    0
end

# ─── CDR3 end in J genes ───

const CDR3_END_REGEX = Dict(
    "IGH" => r"W[GAV]",
    "IGK" => r"FG",
    "IGL" => r"FG",
    "TRA" => r"FG",
    "TRB" => r"FG",
    "TRG" => r"FG",
    "TRD" => r"FG",
)

"""
    cdr3_end_in_j(nt::AbstractString, locus::AbstractString) -> Int

Find CDR3 end position (0-based nucleotide offset) within a J gene sequence.
Returns 0 if not found. Tries all three reading frames.
"""
function cdr3_end_in_j(nt::AbstractString, locus::AbstractString)
    regex = get(CDR3_END_REGEX, locus, nothing)
    regex === nothing && return 0
    for frame in (0, 1, 2)
        frame >= length(nt) && continue
        aa = translate(SubString(nt, frame + 1))
        m = match(regex, aa)
        m !== nothing && return (m.offset - 1) * 3 + frame
    end
    0
end

# ─── Full CDR3 detection in a complete VDJ sequence ───

const CDR3_FIND_IGH = r"[FY][FHVWY]C(?<cdr3>[ADEGIKMNRSTV].{3,31})W[GAV]"
const CDR3_FIND_IGH_ALT = r"C(?<cdr3>.[RK].{3,30})[WF]G.G"
const CDR3_FIND_IGK = r"[FSVY][CFHNVY][CDFGLSW](?<cdr3>.{4,15})[FLV][GRV]"
const CDR3_FIND_IGL = r"[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW](?<cdr3>.{4,15})[FS]G"

"""
    find_cdr3(sequence::AbstractString, locus::AbstractString) -> Tuple{Int,Int}

Find CDR3 region in a nucleotide sequence.
Returns (nt_start, nt_end) tuple (1-based) or (0,0) if not found.
"""
function find_cdr3(sequence::AbstractString, locus::AbstractString)
    regex = if locus == "IGH"
        CDR3_FIND_IGH
    elseif locus == "IGK"
        CDR3_FIND_IGK
    elseif locus == "IGL"
        CDR3_FIND_IGL
    else
        return (0, 0)
    end

    best = (0, 0)
    for offset in (0, 1, 2)
        offset >= length(sequence) && continue
        aa = translate(SubString(sequence, offset + 1))
        m = match(regex, aa)
        if m === nothing && locus == "IGH"
            m = match(CDR3_FIND_IGH_ALT, aa)
        end
        m === nothing && continue
        # Find the 'cdr3' capture group
        idx = findfirst(==("cdr3"), [string(k) for k in keys(m)])
        if idx !== nothing
            start_aa = m.offsets[idx]
            len_aa = length(m.captures[idx])
            nt_start = (start_aa - 1) * 3 + offset + 1
            nt_end = (start_aa - 1 + len_aa) * 3 + offset
            if best == (0, 0) || nt_start < best[1]
                best = (nt_start, nt_end)
            end
        end
    end
    best
end
