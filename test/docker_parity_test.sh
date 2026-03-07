#!/bin/bash
# docker_parity_test.sh — Run Python igdiscover and Julia IgDiscover.jl on the same data,
# then compare outputs for parity.
#
# This script runs inside the Docker container. It:
#   1. Creates a test database from IMGT-style FASTA (with dots in sequences)
#   2. Generates synthetic reads (or uses mounted reads at /data/reads.fasta.gz)
#   3. Runs Python igdiscover
#   4. Runs Julia IgDiscover.jl
#   5. Compares filtered.tsv.gz and new_V_germline outputs
#
# Environment variables:
#   LIMIT           Number of reads to process (default: 0 = all for custom reads, 200 for synthetic)
#   ITERATIONS      Number of iterations (default: 1)
#   SKIP_PYTHON     Set to 1 to skip Python run (reuse existing /results/python_analysis)
#   SKIP_JULIA      Set to 1 to skip Julia run (reuse existing /results/julia_analysis)

set -euo pipefail

LIMIT="${LIMIT:-}"
ITERATIONS="${ITERATIONS:-1}"
SKIP_PYTHON="${SKIP_PYTHON:-0}"
SKIP_JULIA="${SKIP_JULIA:-0}"

RESULTS_DIR="/results"
DB_DIR="${RESULTS_DIR}/database"
PY_DIR="${RESULTS_DIR}/python_analysis"
JL_DIR="${RESULTS_DIR}/julia_analysis"
JL_PROJECT="/igdiscover-jl"

echo "═══════════════════════════════════════════════════"
echo "  IgDiscover Parity Test"
echo "═══════════════════════════════════════════════════"
echo ""
echo "  Python igdiscover: $(igdiscover --version 2>&1 || echo 'not found')"
echo "  Julia:             $(julia --version 2>&1)"
echo "  IgBLAST:           $(igblastn -version 2>&1 | head -1)"
echo "  MUSCLE:            $(muscle -version 2>&1 | head -1 || echo 'v3.x')"
echo ""

# ── Step 1: Determine reads source ──

READS_PATH=""
if [ -f "/data/reads.fasta.gz" ]; then
    READS_PATH="/data/reads.fasta.gz"
    LIMIT="${LIMIT:-0}"
    echo "Using provided reads: $READS_PATH"
elif [ -f "/data/reads.fasta" ]; then
    READS_PATH="/data/reads.fasta"
    LIMIT="${LIMIT:-0}"
    echo "Using provided reads: $READS_PATH"
else
    echo "No reads mounted at /data/reads.fasta.gz"
    echo "Generating synthetic test reads..."
    LIMIT="${LIMIT:-200}"
    julia --project="$JL_PROJECT" "${JL_PROJECT}/test/create_test_data.jl" \
        "${RESULTS_DIR}/synthetic_reads.fasta.gz" \
        "${RESULTS_DIR}/database" \
        "$LIMIT"
    READS_PATH="${RESULTS_DIR}/synthetic_reads.fasta.gz"
    echo "Generated $LIMIT synthetic reads"
fi

echo "  Reads: $READS_PATH"
echo "  Limit: $LIMIT (0 = all)"
echo "  Iterations: $ITERATIONS"
echo ""

# ── Step 2: Create test database (if not already created by synthetic data gen) ──

if [ ! -d "$DB_DIR" ]; then
    echo "── Creating test database ──"
    julia --project="$JL_PROJECT" "${JL_PROJECT}/test/create_test_data.jl" \
        "/dev/null" "$DB_DIR" 0
fi

# Verify database
for gene in V D J; do
    if [ ! -f "${DB_DIR}/${gene}.fasta" ]; then
        echo "ERROR: Missing ${DB_DIR}/${gene}.fasta"
        exit 1
    fi
    COUNT=$(grep -c "^>" "${DB_DIR}/${gene}.fasta" || true)
    echo "  ${gene}.fasta: ${COUNT} sequences"
done

# ── Step 3: Run Python igdiscover ──

if [ "$SKIP_PYTHON" = "1" ] && [ -d "$PY_DIR" ]; then
    echo ""
    echo "── Skipping Python igdiscover (SKIP_PYTHON=1) ──"
else
    echo ""
    echo "── Running Python igdiscover ──"
    rm -rf "$PY_DIR"

    # igdiscover init expects --db with V.fasta, D.fasta, J.fasta
    igdiscover init --db "$DB_DIR" --single-reads "$READS_PATH" "$PY_DIR"

    # Patch config: set limit and iterations
    cd "$PY_DIR"
    if [ "$LIMIT" != "0" ]; then
        sed -i "s/^limit:.*/limit: ${LIMIT}/" igdiscover.yaml
    fi
    sed -i "s/^iterations:.*/iterations: ${ITERATIONS}/" igdiscover.yaml

    echo "  Config:"
    grep -E "^(limit|iterations):" igdiscover.yaml || true
    echo ""

    igdiscover run
    echo "  Python igdiscover completed"
    cd /data
fi

# ── Step 4: Run Julia IgDiscover.jl ──

if [ "$SKIP_JULIA" = "1" ] && [ -d "$JL_DIR" ]; then
    echo ""
    echo "── Skipping Julia IgDiscover.jl (SKIP_JULIA=1) ──"
else
    echo ""
    echo "── Running Julia IgDiscover.jl ──"
    rm -rf "$JL_DIR"

    julia --project="$JL_PROJECT" -e "
        using IgDiscover
        init_analysis(\"${JL_DIR}\", \"${DB_DIR}\", \"${READS_PATH}\")
    "

    # Patch config
    cd "$JL_DIR"
    if [ "$LIMIT" != "0" ]; then
        sed -i "s/^limit = .*/limit = ${LIMIT}/" igdiscover.toml
    fi
    sed -i "s/^iterations = .*/iterations = ${ITERATIONS}/" igdiscover.toml

    echo "  Config:"
    grep -E "^(limit|iterations)" igdiscover.toml || true
    echo ""

    julia --project="$JL_PROJECT" -e "
        using IgDiscover
        run_pipeline(\"${JL_DIR}\")
    "
    echo "  Julia IgDiscover.jl completed"
    cd /data
fi

# ── Step 5: Compare outputs ──

echo ""
echo "── Comparing outputs ──"
julia --project="$JL_PROJECT" "${JL_PROJECT}/test/run_parity_test.jl" "$PY_DIR" "$JL_DIR"
EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "══════════════════════════════════════════════"
    echo "  ✓ PARITY TEST PASSED"
    echo "══════════════════════════════════════════════"
else
    echo "══════════════════════════════════════════════"
    echo "  ✗ PARITY TEST FAILED"
    echo "══════════════════════════════════════════════"
fi

exit $EXIT_CODE
