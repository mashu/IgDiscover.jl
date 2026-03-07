# CDR3 detection — ported from igdiscover/species.py
# Detects CDR3 start within V genes and CDR3 end within J genes.
# Locus is dispatched via the Locus abstract type; string API converts via to_locus.

# ─── Locus dispatch (no string branching in hot path) ───

abstract type Locus end
struct IGH <: Locus end
struct IGK <: Locus end
struct IGL <: Locus end
struct TRA <: Locus end
struct TRB <: Locus end
struct TRG <: Locus end
struct TRD <: Locus end
struct UnknownLocus <: Locus end

const LOCI = Dict{String, Locus}(
    "IGH" => IGH(), "IGK" => IGK(), "IGL" => IGL(),
    "TRA" => TRA(), "TRB" => TRB(), "TRG" => TRG(), "TRD" => TRD(),
)

to_locus(s::AbstractString) = get(LOCI, s, UnknownLocus())

# ─── CDR3 start in V genes ───

const CDR3_START_VH = r"[FY][FHVWY]C([ADEGIKMNRSTV\*]|$)"
const CDR3_START_VH_ALT = r"C(.[RK])"
const CDR3_START_VK = r"[FSVY][CFHNVY][CDFGLSW]"
const CDR3_START_VL = r"[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]"
const CDR3_START_TR = r"[YFH]C"

"""
    cdr3_start_in_v(nt::AbstractString, locus::AbstractString) -> Int

Find CDR3 start position (0-based nucleotide offset) within a V gene sequence.
Returns 0 if not found. Dispatches on Locus type internally.
"""
cdr3_start_in_v(nt::AbstractString, locus::AbstractString) = cdr3_start_in_v(nt, to_locus(locus))

function cdr3_start_in_v(nt::AbstractString, ::IGH)
    aa = translate(nt)
    isempty(aa) && return 0
    cdr3_start_heavy(aa)
end

function cdr3_start_in_v(nt::AbstractString, locus::Union{IGK,IGL,TRG,TRD})
    aa = translate(nt)
    isempty(aa) && return 0
    regex = cdr3_start_regex_v(locus)
    head_len = max(0, length(aa) - 15)
    tail = SubString(aa, head_len + 1)
    m = match(regex, tail)
    m === nothing && return 0
    3 * (head_len + m.offset + length(m.match) - 1)
end
cdr3_start_regex_v(::IGK) = CDR3_START_VK
cdr3_start_regex_v(::IGL) = CDR3_START_VL
cdr3_start_regex_v(::Union{TRG,TRD}) = CDR3_START_TR

function cdr3_start_in_v(nt::AbstractString, ::Union{TRA,TRB})
    aa = translate(nt)
    isempty(aa) && return 0
    head_len = max(0, length(aa) - 8)
    tail = SubString(aa, head_len + 1)
    pos = findfirst('C', tail)
    pos === nothing && return 0
    3 * (head_len + pos)
end

cdr3_start_in_v(::AbstractString, ::UnknownLocus) = 0

function cdr3_start_heavy(aa::AbstractString)
    head_len = max(0, length(aa) - 15)
    tail = SubString(aa, head_len + 1)

    m = match(CDR3_START_VH, tail)
    if m === nothing
        m = match(CDR3_START_VH_ALT, tail)
    end
    m === nothing && return 0

    if length(m.captures) >= 1 && m.captures[1] !== nothing
        cap_offset = m.offsets[1]
        return 3 * (head_len + cap_offset - 1)
    end
    0
end

# ─── CDR3 end in J genes ───

cdr3_end_regex_j(::IGH) = r"W[GAV]"
cdr3_end_regex_j(::IGK) = r"FG"
cdr3_end_regex_j(::IGL) = r"FG"
cdr3_end_regex_j(::TRA) = r"FG"
cdr3_end_regex_j(::TRB) = r"FG"
cdr3_end_regex_j(::TRG) = r"FG"
cdr3_end_regex_j(::TRD) = r"FG"

"""
    cdr3_end_in_j(nt::AbstractString, locus::AbstractString) -> Int

Find CDR3 end position (0-based nucleotide offset) within a J gene sequence.
Returns 0 if not found. Tries all three reading frames. Dispatches on Locus.
"""
cdr3_end_in_j(nt::AbstractString, locus::AbstractString) = cdr3_end_in_j(nt, to_locus(locus))

cdr3_end_in_j(::AbstractString, ::UnknownLocus) = 0

function cdr3_end_in_j(nt::AbstractString, locus::Locus)
    regex = cdr3_end_regex_j(locus)
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
Returns (nt_start, nt_end) tuple (1-based) or (0,0) if not found. Dispatches on Locus.
"""
find_cdr3(sequence::AbstractString, locus::AbstractString) = find_cdr3(sequence, to_locus(locus))

find_cdr3(::AbstractString, ::UnknownLocus) = (0, 0)
find_cdr3(::AbstractString, ::Union{TRA,TRB,TRG,TRD}) = (0, 0)

function find_cdr3(sequence::AbstractString, locus::Union{IGH,IGK,IGL})
    regex_alt = find_cdr3_alt_regex(locus)
    best = (0, 0)
    for offset in (0, 1, 2)
        offset >= length(sequence) && continue
        aa = translate(SubString(sequence, offset + 1))
        m = match(find_cdr3_regex(locus), aa)
        if m === nothing && regex_alt !== nothing
            m = match(regex_alt, aa)
        end
        m === nothing && continue
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

find_cdr3_regex(::IGH) = CDR3_FIND_IGH
find_cdr3_regex(::IGK) = CDR3_FIND_IGK
find_cdr3_regex(::IGL) = CDR3_FIND_IGL
find_cdr3_alt_regex(::IGH) = CDR3_FIND_IGH_ALT
find_cdr3_alt_regex(::Union{IGK,IGL}) = nothing
