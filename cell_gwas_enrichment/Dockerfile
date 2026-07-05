# cell_gwas_enrichment container
#
# Base image: Bioconductor's official docker image, which ships R + BiocManager
# preconfigured with the Bioconductor package repos already registered (this
# avoids most of the pain of installing GenomicRanges/Rsamtools/etc. from
# source). "latest" tracks the current Bioconductor release; pin to a specific
# release tag (e.g. bioconductor/bioconductor_docker:RELEASE_3_20) instead if
# you want a fully reproducible build that won't change out from under you.
FROM bioconductor/bioconductor_docker:latest

LABEL description="Per-cell GWAS variant-effect-magnitude enrichment in ATAC-seq fragments (Signac/Seurat)"

# System libraries needed by transitive dependencies of Seurat/Signac when
# compiling from source (qlcMatrix needs GSL; igraph needs GLPK/GMP; hdf5r
# needs libhdf5). Without these, install.packages() for the affected R
# package can fail silently (it doesn't hard-stop the whole RUN command by
# default) -- that's the most common cause of a "Failed to install: Signac"
# error from the sanity-check steps below despite the install step itself
# reporting success.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgsl-dev \
        libglpk-dev \
        libgmp-dev \
        libhdf5-dev \
        libxml2-dev \
        cmake \
    && rm -rf /var/lib/apt/lists/*

# Bioconductor packages (installed first since they're the slower/trickier part)
RUN Rscript -e 'BiocManager::install(c( \
        "GenomicRanges", \
        "IRanges", \
        "S4Vectors", \
        "GenomeInfoDb", \
        "Rsamtools", \
        "regioneR" \
    ), update = FALSE, ask = FALSE)'
RUN Rscript -e ' \
    required <- c("GenomicRanges","IRanges","S4Vectors","GenomeInfoDb","Rsamtools","regioneR"); \
    missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop("Failed to install: ", paste(missing, collapse = ", "))'

# CRAN packages, split into their own layers so a failure in one doesn't
# force recompiling the others on a rebuild, and so we know exactly which
# package failed instead of finding out only after everything finishes.
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
    install.packages(c("data.table", "ggplot2", "optparse"), repos = "https://cloud.r-project.org")'
RUN Rscript -e ' \
    required <- c("data.table","ggplot2","optparse"); \
    missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop("Failed to install: ", paste(missing, collapse = ", "))'

# Seurat: this is the slow step (large compiled dependency tree -- expect
# 15-20+ minutes on first build).
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
    install.packages("Seurat", repos = "https://cloud.r-project.org")'
RUN Rscript -e 'if (!requireNamespace("Seurat", quietly = TRUE)) stop("Failed to install: Seurat")'

# Signac (depends on Seurat + the Bioconductor packages above).
RUN Rscript -e 'options(Ncpus = parallel::detectCores()); \
    install.packages("Signac", repos = "https://cloud.r-project.org")'
RUN Rscript -e 'if (!requireNamespace("Signac", quietly = TRUE)) stop("Failed to install: Signac")'

WORKDIR /pipeline
COPY R/ ./R/
COPY example_run.R ./
COPY run_pipeline_cli.R ./
COPY tests/ ./tests/
COPY README.md ./
RUN chmod +x run_pipeline_cli.R

# Mount points for user-supplied data (ATAC object, GWAS sumstats) and outputs.
# These are meant to be overridden with -v bind mounts at `docker run` time;
# see README.md, "Running in Docker".
RUN mkdir -p /pipeline/data /pipeline/results

CMD ["bash"]
