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

        @test IgDiscover.hamming_distance("AAAA", "ABCA") == 2

        @test IgDiscover.sequence_hash("ATCG") == IgDiscover.sequence_hash("ATCG")
        @test startswith(IgDiscover.sequence_hash("ATCG"), "S")

        @test IgDiscover.safe_divide(10, 5) == 2.0
        @test IgDiscover.safe_divide(1, 0) == 0.0

        @test IgDiscover.is_same_gene("IGHV1-2*01", "IGHV1-2*02")
        @test !IgDiscover.is_same_gene("IGHV1-2*01", "IGHV1-3*01")
        @test !IgDiscover.is_same_gene("IGHV1", "IGHV2")
    end

    @testset "Configuration" begin
        defaults = joinpath(@__DIR__, "..", "config", "defaults.toml")
        @test isfile(defaults)
        d = IgDiscover.TOML.parsefile(defaults)
        cfg = IgDiscover.parse_config(d)
        @test cfg.iterations == 1
        @test cfg.seed == 1
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

        d = IgDiscover.read_fasta_dict(path)
        @test d["g1"] == "ATCG"
        rm(tmpdir; recursive=true)
    end

    @testset "CDR3 detection" begin
        # Heavy chain CDR3 start: look for [FY][FHVWY]C motif near end
        # Build a fake V gene with the motif at the end
        # YYC is the conserved motif → CDR3 starts after C
        # Amino acid: ... F Y C A ... → CDR3 starts at 'A'
        v_seq = "AAA" ^ 30 * "TTTTATTGT" * "GCT"  # ...FYC then A (CDR3 start)
        pos = IgDiscover.cdr3_start_in_v(v_seq, "IGH")
        @test pos > 0  # should find CDR3 start

        # CDR3 end in J: look for W[GAV] (heavy) or FG (light)
        j_seq = "TGGGCAGGG"  # WA in frame 0
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

    @testset "Clustering" begin
        seqs = ["AAAA", "AAAB", "CCCC", "CCCD"]
        comps = IgDiscover.single_linkage(seqs, (s,t) -> IgDiscover.edit_distance(s,t) <= 1)
        @test length(comps) == 2

        cdr3s = ["AAA", "AAB", "CCC", "CCD", ""]
        j_calls = ["J1", "J1", "J2", "J2", "J1"]
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
            v_call = ["IGHV1", "IGHV2", "", "IGHV3"],
            j_call = ["IGHJ1", "", "IGHJ2", "IGHJ3"],
            stop_codon = ["F", "F", "F", "T"],
            v_support = [1e-5, 1e-1, 1e-5, 1e-5],
            V_covered = [95.0, 85.0, 95.0, 95.0],
            J_covered = [70.0, 70.0, 70.0, 50.0],
            cdr3 = ["AAA", "BBB", "CCC", "DDD"])
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
        ref = IgDiscover.FilterCandidate("ATCG", "g1*01", 10, 5, 3, 100, false, true, true, 4, 1)
        cand = IgDiscover.FilterCandidate("ATCG", "g1*02", 1, 1, 1, 5, false, true, true, 4, 2)

        @test !isempty(IgDiscover.should_discard(IgDiscover.IdenticalSequenceFilter(), ref, cand, false))
        @test !isempty(IgDiscover.should_discard(IgDiscover.ExactRatioFilter(0.5), ref, cand, true))
    end

    @testset "N-tolerant merge" begin
        @test IgDiscover.merge_n_tolerant("ATCG", "ATCG") == "ATCG"
        @test IgDiscover.merge_n_tolerant("ANCG", "ATCG") == "ATCG"
        @test IgDiscover.merge_n_tolerant("ATCG", "GGGC") === nothing
        @test IgDiscover.merge_n_tolerant("ATC", "ATCG") == "ATCG"
    end

    @testset "Reservoir sampling" begin
        data = collect(1:100)
        @test length(IgDiscover.reservoir_sample(data, 10)) == 10
        @test length(IgDiscover.reservoir_sample(data, 200)) == 100
    end

    @testset "Query position mapping" begin
        # Simple case: no gaps, germline pos 3, starts at 1
        pos = IgDiscover.query_position(1, 1, "ATCG", "ATCG", 3)
        @test pos == 3

        # With gap in germline
        pos = IgDiscover.query_position(1, 1, "AT-CG", "ATXCG", 3)
        @test pos >= 3  # query has an extra base at position 3
    end
end
