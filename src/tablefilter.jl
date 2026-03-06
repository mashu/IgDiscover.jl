# Preprocessing filter for assignment tables (produces filtered.tsv.gz)

struct FilterStats
    total::Int
    has_vj_assignment::Int
    has_no_stop::Int
    good_v_evalue::Int
    good_v_coverage::Int
    good_j_coverage::Int
    has_cdr3::Int
end

Base.:(+)(a::FilterStats, b::FilterStats) = FilterStats(
    a.total + b.total,
    a.has_vj_assignment + b.has_vj_assignment,
    a.has_no_stop + b.has_no_stop,
    a.good_v_evalue + b.good_v_evalue,
    a.good_v_coverage + b.good_v_coverage,
    a.good_j_coverage + b.good_j_coverage,
    a.has_cdr3 + b.has_cdr3)

"""
    filter_table(df, pf) -> (DataFrame, FilterStats)

Apply preprocessing filter: require V+J assignment, no stop codons, coverage/evalue thresholds.
"""
function filter_table(df::DataFrame, pf::PreprocessingFilter)
    total = nrow(df)

    filtered = df[df.v_call .!= "" .&& df.j_call .!= "", :]
    has_vj = nrow(filtered)

    if hasproperty(filtered, :stop_codon)
        no_stop = [let s = lowercase(string(coalesce(x, "")))
            s == "f" || s == "false"
        end for x in filtered.stop_codon]
        filtered = filtered[no_stop, :]
    end
    has_no_stop = nrow(filtered)

    if hasproperty(filtered, :v_support)
        filtered = filtered[coalesce.(filtered.v_support, Inf) .<= pf.v_evalue, :]
    end
    good_v_eval = nrow(filtered)

    if hasproperty(filtered, :V_covered)
        filtered = filtered[coalesce.(filtered.V_covered, 0.0) .>= pf.v_coverage, :]
    end
    good_v_cov = nrow(filtered)

    if hasproperty(filtered, :J_covered)
        filtered = filtered[coalesce.(filtered.J_covered, 0.0) .>= pf.j_coverage, :]
    end
    good_j_cov = nrow(filtered)

    has_cdr3 = hasproperty(filtered, :cdr3) ? count(!isempty, filtered.cdr3) : 0

    stats = FilterStats(total, has_vj, has_no_stop, good_v_eval, good_v_cov, good_j_cov, has_cdr3)
    @info "Filter: $total → $has_vj (V+J) → $has_no_stop (no stop) → $good_v_eval (evalue) → $good_v_cov (V cov) → $good_j_cov (J cov)"

    (filtered, stats)
end

"""
    filter_table(input_path, output_path, pf; stats_path) -> (DataFrame, FilterStats)

File-based variant: reads input, filters, writes output and optional JSON stats.
"""
function filter_table(
    input_path::AbstractString,
    output_path::AbstractString,
    pf::PreprocessingFilter;
    stats_path::String="",
)
    df = read_assignments(input_path)
    filtered, stats = filter_table(df, pf)

    endswith(output_path, ".gz") ? write_table_gz(output_path, filtered) :
                                   write_table(output_path, filtered)

    if !isempty(stats_path)
        open(stats_path, "w") do io
            JSON3.write(io, Dict(
                "total"             => stats.total,
                "has_vj_assignment" => stats.has_vj_assignment,
                "has_no_stop"       => stats.has_no_stop,
                "good_v_evalue"     => stats.good_v_evalue,
                "good_v_coverage"   => stats.good_v_coverage,
                "good_j_coverage"   => stats.good_j_coverage,
                "has_cdr3"          => stats.has_cdr3,
            ))
        end
    end
    (filtered, stats)
end
