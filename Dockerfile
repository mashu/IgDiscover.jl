# IgDiscover.jl Parity Test Environment
#
# Installs both Python igdiscover (via conda/bioconda) and Julia IgDiscover.jl,
# then provides infrastructure to run both on the same data and compare outputs.
#
# Build:
#   docker build -t igdiscover-parity .
#
# Run parity test (uses built-in synthetic test data):
#   docker run --rm igdiscover-parity
#
# Run with custom reads:
#   docker run --rm -v /path/to/reads.fasta.gz:/data/reads.fasta.gz igdiscover-parity
#
# Interactive:
#   docker run --rm -it igdiscover-parity bash

FROM mambaorg/micromamba:1.5-jammy AS base

USER root

# ── System dependencies ──
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ca-certificates git curl \
    && rm -rf /var/lib/apt/lists/*

# ── Python igdiscover via conda (bioconda bundles igblast, muscle, pear) ──
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        igdiscover=0.15.1 \
        python=3.10 \
    && micromamba clean --all --yes

# Make conda env available
ENV PATH="/opt/conda/bin:$PATH"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# Verify Python igdiscover works
RUN igdiscover --version

# ── Julia installation ──
ARG JULIA_VERSION=1.11.3
RUN JULIA_MAJOR=$(echo $JULIA_VERSION | cut -d. -f1-2) && \
    wget -q "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_MAJOR}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" && \
    tar xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C /opt && \
    ln -s /opt/julia-${JULIA_VERSION}/bin/julia /usr/local/bin/julia && \
    rm julia-${JULIA_VERSION}-linux-x86_64.tar.gz

RUN julia --version

# ── IgDiscover.jl project setup ──
WORKDIR /igdiscover-jl
COPY Project.toml .
COPY src/ src/
COPY config/ config/
COPY test/ test/
COPY bin/ bin/

# Precompile Julia dependencies (slow, but only done once during build)
RUN julia --project=. -e ' \
    using Pkg; \
    Pkg.instantiate(); \
    Pkg.precompile(); \
    '

# ── Test infrastructure ──
COPY test/docker_parity_test.sh /usr/local/bin/run_parity_test
RUN chmod +x /usr/local/bin/run_parity_test

# Working directory for test runs
RUN mkdir -p /data /results
WORKDIR /data

# Default: run the parity test
CMD ["run_parity_test"]
