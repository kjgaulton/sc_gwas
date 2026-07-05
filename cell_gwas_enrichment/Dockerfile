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

# Bioconductor packages (installed first since they're the slower/trickier part)
RUN Rscript -e 'BiocManager::install(c( \
        "GenomicRanges", \
        "IRanges", \
        "S4Vectors", \
        "GenomeInfoDb", \
        "Rsamtools", \
        "regioneR" \
    ), update = FALSE, ask = FALSE)'

# CRAN packages, including Seurat + Signac. These pull in a large number of
# transitive dependencies (Matrix, RcppArmadillo, irlba, uwot, cowplot, ...)
# and will take a while to build from source on first image build.
RUN Rscript -e 'install.packages(c( \
        "data.table", \
        "ggplot2", \
        "optparse", \
        "Seurat", \
        "Signac" \
    ), repos = "https://cloud.r-project.org")'

# Sanity check: fail the build early (with a clear error) if any package
# didn't actually install, rather than surfacing a confusing error later
# at pipeline run time.
RUN Rscript -e ' \
    required <- c("GenomicRanges","IRanges","S4Vectors","GenomeInfoDb","Rsamtools", \
                  "regioneR","data.table","ggplot2","optparse","Seurat","Signac"); \
    missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop("Failed to install: ", paste(missing, collapse = ", "))'

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
