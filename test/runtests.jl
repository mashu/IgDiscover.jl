using Test
using IgDiscover
using DataFrames

@testset "IgDiscover.jl" begin

    @testset "DNA utilities" begin
        @test IgDiscover.translate("ATGAAAGGG") == "MKG"
        @test IgDiscover.translate("TAATAGTGA") == "***"
        @test IgDiscover.reverse_complement("ATCG") == "CGAT"
        @test IgDiscover.has_stop("TAA")
        @test !IgDiscover.has_stop("GGG")
        @test !IgDiscover.has_stop("GGGAC")

        @test IgDiscover.edit_distance("kitten", "sitting") == 3
        @test IgDiscover.edit_distance("ABC", "ABC") == 0
        @test IgDiscover.edit_distance("ABC", "ABD") == 1
        @test IgDiscover.edit_distance("ABC", "XYZ"; maxdiff=1) == 2
        @test IgDiscover.edit_distance("", "") == 0
        @test IgDiscover.edit_distance("ABC", "") == 3
        @test IgDiscover.edit_distance("", "ABC") == 3

        @test IgDiscover.hamming_distance("AAAA", "ABCA") == 2

        @test IgDiscover.sequence_hash("ATCG") == IgDiscover.sequence_hash("ATCG")
        @test startswith(IgDiscover.sequence_hash("ATCG"), "S")

        @test IgDiscover.safe_divide(10, 5) == 2.0
        @test IgDiscover.safe_divide(1, 0) == 0.0

        @test IgDiscover.is_same_gene("IGHV1-2*01", "IGHV1-2*02")
        @test !IgDiscover.is_same_gene("IGHV1-2*01", "IGHV1-3*01")
        @test !IgDiscover.is_same_gene("IGHV1", "IGHV2")
    end

    @testset "Edit distance thread safety" begin
        results = Vector{Int}(undef, 100)
        Threads.@threads for i in 1:100
            results[i] = IgDiscover.edit_distance("ABCDEF", "XBCDEY")
        end
        @test all(==(2), results)

        results2 = Vector{Int}(undef, 50)
        Threads.@threads for i in 1:50
            s = repeat("A", i)
            t = repeat("A", i)
            results2[i] = IgDiscover.edit_distance(s, t)
        end
        @test all(==(0), results2)
    end

    @testset "Edit distance with maxdiff" begin
        @test IgDiscover.edit_distance("AAAA", "ZZZZ"; maxdiff=2) == 3
        @test IgDiscover.edit_distance("AAAA", "AABA"; maxdiff=5) == 1
        @test IgDiscover.edit_distance("A", "ABCD"; maxdiff=1) == 2
    end

    @testset "Configuration" begin
        defaults = joinpath(@__DIR__, "..", "config", "defaults.toml")
        @test isfile(defaults)
        d = IgDiscover.TOML.parsefile(defaults)
        cfg = IgDiscover.parse_config(d)
        @test cfg.iterations == 1
        @test cfg.seed == 1
        @test cfg.consensus_threshold == 60.0
        @test cfg.preprocessing_filter.v_coverage == 90.0
        @test cfg.germline_filter.unique_cdr3s == 5
        @test cfg.pre_germline_filter.unique_cdr3s == 2
        @test cfg.j_discovery.propagate
    end

    @testset "FASTA I/O" begin
        tmpdir = mktempdir()
        path = joinpath(tmpdir, "test.fasta")
        records = [IgDiscover.FastaRecord("g1", "ATCG"),
                   IgDiscover.FastaRecord("g2", "GGCC")]
        IgDiscover.write_fasta(path, records)
        loaded = IgDiscover.read_fasta(path)
        @test length(loaded) == 2
        @test loaded[1].name == "g1"
        @test loaded[1].sequence == "ATCG"

        gz_path = joinpath(tmpdir, "test.fasta.gz")
        IgDiscover.write_fasta_gz(gz_path, records)
        loaded_gz = IgDiscover.read_fasta(gz_path)
        @test length(loaded_gz) == 2
        @test loaded_gz[1].name == "g1"

        d = IgDiscover.read_fasta_dict(path)
        @test d["g1"] == "ATCG"

        loaded_lim = IgDiscover.read_fasta(path; limit=1)
        @test length(loaded_lim) == 1

        rm(tmpdir; recursive=true)
    end

    @testset "IMGT sanitization" begin
        # allele_name_from_header
        @test IgDiscover.allele_name_from_header("J00256|IGHJ1*01|Homo sapiens|F|") == "IGHJ1*01"
        @test IgDiscover.allele_name_from_header("M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|") == "IGHV1-18*01"
        @test IgDiscover.allele_name_from_header("IGHV1-18*01") == "IGHV1-18*01"
        @test IgDiscover.allele_name_from_header("simple") == "simple"

        # sanitize_imgt_sequence removes dots
        @test IgDiscover.sanitize_imgt_sequence("ATG...CCC.GGG") == "ATGCCCGGG"
        @test IgDiscover.sanitize_imgt_sequence("ATCG") == "ATCG"
        @test IgDiscover.sanitize_imgt_sequence("...") == ""
        @test IgDiscover.sanitize_imgt_sequence("atg..ccc") == "ATGCCC"

        # sanitize_imgt_record
        rec = IgDiscover.FastaRecord("M99641|IGHV1-18*01|Homo sapiens|F|", "atg...ccc.ggg")
        cleaned = IgDiscover.sanitize_imgt_record(rec)
        @test cleaned.name == "IGHV1-18*01"
        @test cleaned.sequence == "ATGCCCGGG"

        # Round-trip through write_sanitized_imgt
        tmpdir = mktempdir()
        imgt_path = joinpath(tmpdir, "imgt.fasta")
        out_path = joinpath(tmpdir, "clean.fasta")

        open(imgt_path, "w") do io
            println(io, ">M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |")
            println(io, "cagg...ttcagctggtgcag...tctggagct...gaggtg...aagaagcctggggcc")
            println(io, ">X62106|IGHV1-18*04|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |")
            println(io, "cagg.ttcagctggtgcag.tctggagct")
        end

        IgDiscover.write_sanitized_imgt(imgt_path, out_path)
        result = IgDiscover.read_fasta(out_path)
        @test length(result) == 2
        @test result[1].name == "IGHV1-18*01"
        @test !occursin('.', result[1].sequence)
        @test result[1].sequence == "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCC"
        @test result[2].name == "IGHV1-18*04"
        @test result[2].sequence == "CAGGTTCAGCTGGTGCAGTCTGGAGCT"

        rm(tmpdir; recursive=true)
    end

    @testset "CDR3 detection" begin
        v_seq = "AAA" ^ 30 * "TTTTATTGT" * "GCT"  # ...FYC then A (CDR3 start)
        pos = IgDiscover.cdr3_start_in_v(v_seq, "IGH")
        @test pos > 0

        j_seq = "TGGGCAGGG"  # WA motif
        pos_j = IgDiscover.cdr3_end_in_j(j_seq, "IGH")
        @test pos_j >= 0
    end

    @testset "Alignment" begin
        aln = IgDiscover.align_affine("AAA", "AGA")
        @test aln.errors == 1
        @test IgDiscover.describe_nt_change("AAA", "AGA") == "2A>G"
        @test IgDiscover.describe_nt_change("AAGG", "AATTGG") == "2_3insTT"
        @test IgDiscover.describe_nt_change("AATTGG", "AAGG") == "3_4delTT"
    end

    @testset "NucleotideCounter" begin
        nc = IgDiscover.NucleotideCounter()
        IgDiscover.count!(nc, 'A')
        IgDiscover.count!(nc, 'A')
        IgDiscover.count!(nc, 'C')
        IgDiscover.count!(nc, '-')
        ch, freq = IgDiscover.best_base(nc)
        @test ch == 'A'
        @test freq == 2

        IgDiscover.reset!(nc)
        @test nc.a == 0
        @test nc.gap == 0
    end

    @testset "Consensus sequence" begin
        seqs = ["AAAA", "AABA", "AACA"]
        cons = IgDiscover.consensus_sequence(seqs; threshold=0.5)
        @test length(cons) == 4
        @test cons[1] == 'A'

        @test IgDiscover.consensus_sequence(String[]) == ""
        @test IgDiscover.consensus_sequence(["ATCG"]; threshold=0.5) == "ATCG"

        seqs2 = ["ATCG", "ATCG", "ATCG"]
        @test IgDiscover.consensus_sequence(seqs2; threshold=0.7) == "ATCG"
    end

    @testset "Clustering" begin
        seqs = ["AAAA", "AAAB", "CCCC", "CCCD"]
        comps = IgDiscover.single_linkage(seqs, (s, t) -> IgDiscover.edit_distance(s, t) <= 1)
        @test length(comps) == 2

        cdr3s   = ["AAA", "AAB", "CCC", "CCD", ""]
        j_calls = ["J1",  "J1",  "J2",  "J2",  "J1"]
        n = IgDiscover.count_clonotypes(cdr3s, j_calls; max_distance=1)
        @test n == 2
    end

    @testset "Header parsing" begin
        name, count, bc = IgDiscover.parse_header("myread;size=12;barcode=ACG;")
        @test name == "myread"
        @test count == 12
        @test bc == "ACG"

        name2, count2, bc2 = IgDiscover.parse_header("simple_name")
        @test name2 == "simple_name"
        @test count2 == 0
        @test bc2 == ""
    end

    @testset "Table filtering" begin
        df = DataFrame(
            v_call     = ["IGHV1", "IGHV2", "", "IGHV3"],
            j_call     = ["IGHJ1", "", "IGHJ2", "IGHJ3"],
            stop_codon = ["F", "F", "F", "T"],
            v_support  = [1e-5, 1e-1, 1e-5, 1e-5],
            V_covered  = [95.0, 85.0, 95.0, 95.0],
            J_covered  = [70.0, 70.0, 70.0, 50.0],
            cdr3       = ["AAA", "BBB", "CCC", "DDD"])
        pf = IgDiscover.PreprocessingFilter(90.0, 60.0, 1e-3)
        filtered, stats = IgDiscover.filter_table(df, pf)
        @test stats.total == 4
        @test stats.has_vj_assignment == 2
        @test nrow(filtered) == 1
    end

    @testset "NameGenerator" begin
        ng = IgDiscover.NameGenerator()
        @test ng("g1") == "g1"
        @test ng("g1") == "g1A"
        @test ng("g1") == "g1B"
        @test ng("g2") == "g2"
    end

    @testset "Germline filter dispatch" begin
        ref  = IgDiscover.FilterCandidate("ATCG", "g1*01", 10, 5, 3, 100, false, true, true, 4, 1)
        cand = IgDiscover.FilterCandidate("ATCG", "g1*02",  1, 1, 1,   5, false, true, true, 4, 2)

        @test !isempty(IgDiscover.should_discard(IgDiscover.IdenticalSequenceFilter(), ref, cand, false))
        @test !isempty(IgDiscover.should_discard(IgDiscover.ExactRatioFilter(0.5), ref, cand, true))

        cand_wl = IgDiscover.FilterCandidate("ATCG", "g1*02", 1, 1, 1, 5, true, true, true, 4, 2)
        @test isempty(IgDiscover.should_discard(IgDiscover.IdenticalSequenceFilter(), ref, cand_wl, false))
    end

    @testset "N-tolerant merge" begin
        @test IgDiscover.merge_n_tolerant("ATCG", "ATCG") == "ATCG"
        @test IgDiscover.merge_n_tolerant("ANCG", "ATCG") == "ATCG"
        @test IgDiscover.merge_n_tolerant("ATCG", "GGGC") === nothing
        @test IgDiscover.merge_n_tolerant("ATC",  "ATCG") == "ATCG"
    end

    @testset "Reservoir sampling" begin
        data = collect(1:100)
        @test length(IgDiscover.reservoir_sample(data, 10))  == 10
        @test length(IgDiscover.reservoir_sample(data, 200)) == 100
    end

    @testset "Query position mapping" begin
        pos = IgDiscover.query_position(1, 1, "ATCG", "ATCG", 3)
        @test pos == 3

        pos = IgDiscover.query_position(1, 1, "AT-CG", "ATXCG", 3)
        @test pos >= 3
    end

    @testset "Barcode extraction" begin
        rec = IgDiscover.FastaRecord("r1", "AAAAATCGATCG")
        bc, unbarcoded = IgDiscover.extract_barcode(rec, 5)
        @test bc == "AAAAA"
        @test unbarcoded.sequence == "TCGATCG"

        bc2, unbarcoded2 = IgDiscover.extract_barcode(rec, -4)
        @test bc2 == "ATCG"
        @test unbarcoded2.sequence == "AAAAATCG"
    end

    @testset "Pseudo-CDR3" begin
        seq = "ATCGATCGATCGATCGATCG"  # length 20
        p = IgDiscover.pseudo_cdr3(seq, 15, 5)
        @test !isempty(p)
        @test length(p) == 10
    end

    @testset "Trim leading G" begin
        @test IgDiscover.trim_leading_g("GGGAATCG") == "AATCG"
        @test IgDiscover.trim_leading_g("AATCG")    == "AATCG"
        @test IgDiscover.trim_leading_g("GGGG")      == ""
    end

    @testset "J discovery filter helpers" begin
        j1 = IgDiscover.JCandidate("J1*01", "J1*01", 100, 10, 0.0, 0.0, "ATCGATCG")
        j2 = IgDiscover.JCandidate("J1*02", "J1*02",  10,  5, 0.0, 0.0, "ATCGATCC")
        j3 = IgDiscover.JCandidate("J1*03", "J1*03",   1,  1, 0.0, 0.0, "ATCGATCA")

        filtered = IgDiscover.filter_j_alleles([j1, j2, j3], 0.2)
        @test j1 in filtered
        @test !(j3 in filtered)
    end

    @testset "Tallies" begin
        t = IgDiscover.tallies(["a", "b", "a", "c", "a"])
        @test t["a"] == 3
        @test t["b"] == 1
        @test t["c"] == 1
    end

    @testset "Rename genes" begin
        tmpdir = mktempdir()
        db_path = joinpath(tmpdir, "db.fasta")
        disc_path = joinpath(tmpdir, "disc.fasta")
        out_path = joinpath(tmpdir, "out.fasta")

        IgDiscover.write_fasta(db_path, [
            IgDiscover.FastaRecord("IGHV1-2*01", "ATCGATCG"),
            IgDiscover.FastaRecord("IGHV1-3*01", "GGCCGGCC"),
        ])
        IgDiscover.write_fasta(disc_path, [
            IgDiscover.FastaRecord("IGHV1-2_S1234", "ATCGATCG"),
            IgDiscover.FastaRecord("IGHV1-2_S5678", "ATCGATCA"),
        ])

        IgDiscover.rename_genes(disc_path, db_path, out_path)
        result = IgDiscover.read_fasta(out_path)
        @test length(result) == 2
        @test any(r -> startswith(r.name, "IGHV1-2*"), result)
        rm(tmpdir; recursive=true)
    end

    @testset "Distance matrix" begin
        seqs = ["AAAA", "AABA", "AAAA"]
        M, clusters = IgDiscover.cluster_sequences(seqs; minsize=1)
        @test size(M) == (3, 3)
        @test M[1, 1] == 0.0
        @test M[1, 3] == 0.0
        @test M[1, 2] > 0.0
    end

    @testset "Validate FASTA" begin
        tmpdir = mktempdir()

        path = joinpath(tmpdir, "valid.fasta")
        IgDiscover.write_fasta(path, [
            IgDiscover.FastaRecord("a", "ATCG"),
            IgDiscover.FastaRecord("b", "GGCC"),
        ])
        IgDiscover.validate_fasta(path)

        dup_path = joinpath(tmpdir, "dup.fasta")
        IgDiscover.write_fasta(dup_path, [
            IgDiscover.FastaRecord("a", "ATCG"),
            IgDiscover.FastaRecord("a", "GGCC"),
        ])
        @test_throws ErrorException IgDiscover.validate_fasta(dup_path)

        rm(tmpdir; recursive=true)
    end

end
