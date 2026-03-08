# build/snoop_workload.jl
#
# Exercises IgDiscover.jl code paths for PackageCompiler to trace.
# This runs during sysimage creation and ensures all critical methods
# are compiled into the image.

using IgDiscover
using DataFrames

tmpdir = mktempdir()

# ── FASTA I/O ──

fasta_path = joinpath(tmpdir, "test.fasta")
open(fasta_path, "w") do io
    for i in 1:20
        seq = join(rand("ACGT", 300))
        println(io, ">seq$i")
        println(io, seq)
    end
end

records = read_fasta(fasta_path)
db = read_fasta_dict(fasta_path)

out_path = joinpath(tmpdir, "out.fasta")
write_fasta(out_path, records)

gz_path = joinpath(tmpdir, "out.fasta.gz")
write_fasta_gz(gz_path, records)
read_fasta(gz_path)

# ── DNA utilities ──

for (s, t) in [("ATCGATCG", "ATCGATCA"), ("AAAA", "BBBB"), ("", "ABC")]
    IgDiscover.edit_distance(s, t)
    length(s) > 0 && length(t) > 0 && IgDiscover.edit_distance(s, t; maxdiff=2)
end
IgDiscover.hamming_distance("ATCG", "ATCA")
IgDiscover.translate("ATGATGATGATG")
IgDiscover.reverse_complement("ATCGATCG")
IgDiscover.has_stop("ATGATGTAA")
IgDiscover.sequence_hash("ATCGATCG")
IgDiscover.unique_name("IGHV1-2*01", "ATCGATCG")

# ── CDR3 detection ──

v_seq = "AAA" ^ 30 * "TTTTATTGTGCT"
IgDiscover.cdr3_start_in_v(v_seq, "IGH")
IgDiscover.cdr3_start_in_v(v_seq, "IGK")
IgDiscover.cdr3_end_in_j("TGGGCAGGG", "IGH")

# ── Config ──

defaults_path = joinpath(dirname(dirname(pathof(IgDiscover))), "config", "defaults.toml")
if isfile(defaults_path)
    raw = IgDiscover.TOML.parsefile(defaults_path)
    IgDiscover.parse_config(raw)
end

# ── Consensus ──

seqs = [join(rand("ACGT", 100)) for _ in 1:10]
seqs_similar = [seqs[1] for _ in 1:8]
push!(seqs_similar, replace(seqs[1], r"^.{3}" => "NNN"))
IgDiscover.consensus_sequence(seqs_similar; threshold=0.6)

# ── Alignment ──

IgDiscover.align_affine("ATCGATCG", "ATCAATCG")
IgDiscover.describe_nt_change("ATCGATCG", "ATCAATCG")

# ── Clustering ──

test_seqs = ["AAAA", "AABA", "CCCC", "CCCD"]
IgDiscover.single_linkage(test_seqs,
    (s, t) -> IgDiscover.edit_distance(s, t) <= 1)
IgDiscover.count_clonotypes(
    ["AAA", "AAB", "CCC"], ["J1", "J1", "J2"]; max_distance=1)
IgDiscover.tallies(["a", "b", "a", "c"])

# ── Header parsing ──

IgDiscover.parse_header("read1;size=5;barcode=ACG;")

# ── IMGT sanitization ──

IgDiscover.sanitize_imgt_sequence("ATG...CCC.GGG")
IgDiscover.allele_name_from_header("SYN001|IGHV1-18*01|Synthetic")

# ── Table operations ──

df = DataFrame(
    v_call = ["IGHV1", "IGHV2"],
    j_call = ["IGHJ1", "IGHJ2"],
    cdr3 = ["AATGT", "CCCGT"],
    V_SHM = [1.0, 2.0],
    stop_codon = ["F", "F"],
    v_support = [1e-5, 1e-4],
    V_covered = [95.0, 90.0],
    J_covered = [70.0, 65.0],
)
pf = IgDiscover.PreprocessingFilter(90.0, 60.0, 1e-3)
IgDiscover.filter_table(df, pf)

# ── Germline filter ──

IgDiscover.FilterCandidate("ATCG", "g1*01", 10, 5, 3, 100, false, true, true, 4, 1)
IgDiscover.IdenticalSequenceFilter()
IgDiscover.CrossMappingFilter(0.02)

# ── Clonotypes ──

ct_df = DataFrame(
    v_call = ["IGHV1", "IGHV1", "IGHV2"],
    j_call = ["IGHJ1", "IGHJ1", "IGHJ2"],
    cdr3   = ["AATGT", "AATGT", "CCCGT"],
    V_SHM  = [1.0, 2.0, 0.5],
)
IgDiscover.call_clonotypes(ct_df)

# ── Cleanup ──

rm(tmpdir; recursive=true, force=true)

# ── App entry point (for standalone build) ──

let orig = copy(ARGS)
    empty!(ARGS)
    push!(ARGS, "version")
    IgDiscover.julia_main()
    empty!(ARGS)
    append!(ARGS, orig)
end

@info "Snoop workload completed successfully"
