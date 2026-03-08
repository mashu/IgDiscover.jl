using Test
using DataFrames
using IgDiscover

@testset "CDR3CoreSlice and is_similar_with_junction" begin
    @test is_similar_with_junction("ATCGATCG", "ATCGATCG", 1.0, nothing) == true
    @test is_similar_with_junction("ATCGATCG", "ATCGATCA", 1.0, nothing) == true
    @test is_similar_with_junction("ATCGATCG", "ATCGATCA", 0.0, nothing) == false
    @test is_similar_with_junction("ATCGATCG", "ATCG", 1.0, nothing) == false

    # Fractional mismatches
    @test is_similar_with_junction("ATCGATCG", "ATCGATCA", 0.15, nothing) == true
    @test is_similar_with_junction("ATCGATCG", "AACAATCA", 0.1, nothing) == false

    # With core slice
    core = CDR3CoreSlice(2, 6)
    @test is_similar_with_junction("ABCDEFGH", "ABCXEFGH", 2.0, core) == true
    @test is_similar_with_junction("ABCDEFGH", "XBCDEFXH", 2.0, core) == false
end

@testset "ClonoQuery" begin
    ref = DataFrame(
        sequence_id = ["r1", "r2", "r3", "r4"],
        v_call = ["IGHV1", "IGHV1", "IGHV2", "IGHV1"],
        j_call = ["IGHJ1", "IGHJ1", "IGHJ1", "IGHJ2"],
        cdr3 = ["ATCGATCG", "ATCGATCA", "GGCCGGCC", "ATCGATCG"],
        count = [10, 5, 3, 1],
    )
    query = DataFrame(
        sequence_id = ["q1", "q2"],
        v_call = ["IGHV1", "IGHV2"],
        j_call = ["IGHJ1", "IGHJ1"],
        cdr3 = ["ATCGATCG", "TTTTTTTT"],
    )

    cq = ClonoQuery(mismatches=1.0)
    results = cq(query, ref)
    @test length(results) == 2
    @test results[1].query_id == "q1"
    @test results[1].match_count == 2
    @test results[2].match_count == 0

    cq2 = ClonoQuery(mismatches=1.0, min_count=6)
    results2 = cq2(query, ref)
    @test results2[1].match_count == 1
end

@testset "AlleleUsage" begin
    table = DataFrame(
        v_call = ["IGHV1*01", "IGHV1*01", "IGHV2*01", "IGHV1*01", "IGHV2*01"],
        j_call = ["IGHJ1*01", "IGHJ2*01", "IGHJ1*01", "IGHJ1*01", "IGHJ2*01"],
        V_errors = [0, 0, 0, 0, 0],
        J_errors = [0, 0, 0, 0, 0],
    )
    au = AlleleUsage(x_gene=:V, y_gene=:J)
    matrix = au(table, ["IGHJ1*01", "IGHJ2*01"])
    @test nrow(matrix) == 2
    v1_row = findfirst(==("IGHV1*01"), matrix.gene)
    @test matrix[v1_row, Symbol("IGHJ1*01")] == 2
    @test matrix[v1_row, Symbol("IGHJ2*01")] == 1
end

@testset "ExpressionCounter" begin
    table = DataFrame(
        v_call = ["IGHV1*01", "IGHV1*01", "IGHV1*02", "IGHV2*01"],
        j_call = ["IGHJ1", "IGHJ2", "IGHJ1", "IGHJ1"],
        d_call = ["IGHD1", "IGHD1", "IGHD2", "IGHD1"],
    )
    counter = ExpressionCounter(gene=:V)
    counts = counter(table)
    @test nrow(counts) == 3
    @test counts.count[findfirst(==("IGHV1*01"), counts.gene)] == 2

    counter_ar = ExpressionCounter(gene=:V, allele_ratio=0.6)
    counts_ar = counter_ar(table)
    @test !("IGHV1*02" in counts_ar.gene)

    counts2 = counter(table; gene_names=["IGHV1*01", "IGHV1*02", "IGHV3*01"])
    @test nrow(counts2) == 3
    @test counts2.count[findfirst(==("IGHV3*01"), counts2.gene)] == 0
end

@testset "natural_sort_key" begin
    names = ["IGHV1-10*01", "IGHV1-2*01", "IGHV1-1*01"]
    @test sort(names; by=natural_sort_key) == ["IGHV1-1*01", "IGHV1-2*01", "IGHV1-10*01"]
end

@testset "filter_by_allele_ratio" begin
    df = DataFrame(gene=["A*01", "A*02", "B*01"], count=[100, 40, 50])
    result = filter_by_allele_ratio(df, 0.5)
    @test nrow(result) == 2
    @test "A*02" ∉ result.gene  # 40/100 = 0.4 < 0.5
end

@testset "DatabaseComparator" begin
    a = [FastaRecord("g1", "ATCGATCG"), FastaRecord("g2", "GGCCGGCC"), FastaRecord("g3", "TTTTTTTT")]
    b = [FastaRecord("g1", "ATCGATCG"), FastaRecord("g2x", "GGCCGGCC"), FastaRecord("g4", "AAAAAAAA")]

    diff = DatabaseComparator()(a, b)
    @test length(diff.identical) == 2
    @test length(diff.renamed) == 1
    @test length(diff.only_a) == 1
    @test length(diff.only_b) == 1
    @test exit_code(diff) == 1

    diff2 = DatabaseComparator()(a[1:2], a[1:2])
    @test exit_code(diff2) == 0
end

@testset "check_duplicate_names" begin
    recs = [FastaRecord("a", "ATCG"), FastaRecord("b", "GGCC"), FastaRecord("a", "TTTT")]
    @test check_duplicate_names(recs) == ["a"]
end

@testset "check_duplicate_sequences" begin
    recs = [FastaRecord("a", "ATCG"), FastaRecord("b", "ATCG")]
    dups = check_duplicate_sequences(recs)
    @test length(dups) == 1
    @test dups[1] == ("b", "a")
end

@testset "UpstreamParams construction" begin
    p = UpstreamParams()
    @test p.max_v_errors == 1.0
    @test p.part == :UTR_leader
    @test_throws ErrorException UpstreamParams(part=:invalid)
end

@testset "format_diff output" begin
    diff = DatabaseComparator()([FastaRecord("g1", "ATCG")], [FastaRecord("g1", "ATCG")])
    buf = IOBuffer()
    format_diff(buf, diff)
    output = String(take!(buf))
    @test occursin("0 lost", output)
    @test occursin("0 gained", output)
end

@testset "HaplotypeEntry and HaplotypePair" begin
    entry = HaplotypeEntry("IGHV1*01", "IGHV1*02", :heterozygous, 100, 80)
    @test entry.classification == :heterozygous

    hp = HaplotypePair([entry], 'V', "IGHJ1*01", "IGHJ1*02")
    tsv = format_tsv(hp)
    @test occursin("haplotype1", tsv)
    @test occursin("IGHV1*01", tsv)

    switch_haplotype!(hp)
    @test hp.entries[1].hap1 == "IGHV1*02"
    @test hp.entries[1].hap2 == "IGHV1*01"
end

@testset "ClusterPlotter construction" begin
    plotter = ClusterPlotter(min_group_size=50)
    @test plotter.min_group_size == 50
    @test plotter.gene_type == :V
end

@testset "subsample" begin
    pop = ["a", "b", "c", "d", "e", "f", "g"]
    sample = IgDiscover.subsample(pop, 3)
    @test length(sample) == 3
    @test all(s -> s in pop, sample)

    sample2 = IgDiscover.subsample(pop, 100)
    @test length(sample2) == length(pop)
end

@testset "pairwise_distance_matrix" begin
    seqs = ["ATCG", "ATCA", "GGCC"]
    mat = IgDiscover.pairwise_distance_matrix(seqs)
    @test size(mat) == (3, 3)
    @test mat[1, 1] == 0.0
    @test mat[1, 2] == mat[2, 1]
    @test mat[1, 2] > 0.0
end

@testset "render_heatmap" begin
    mat = [0.0 1.0; 1.0 0.0]
    result = ClusterResult("TestGene", 2, 1, [1, 1], mat)
    output = render_heatmap(result)
    @test occursin("TestGene", output)
    @test occursin("2 sequences", output)
end

@testset "apply_filters dispatch" begin
    table = DataFrame(
        v_call = ["A", "B"],
        V_errors = [0, 1],
        J_errors = [0, 0],
    )
    au = AlleleUsage(x_gene=:V, y_gene=:J)
    filtered = apply_filters(table, au)
    @test nrow(filtered) == 1

    ha = HaplotypeAnalyzer()
    filtered2 = apply_filters(table, ha)
    @test nrow(filtered2) == 1
end
