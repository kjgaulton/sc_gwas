#!/usr/bin/env Rscript
#' run_pipeline_cli.R
#'
#' Command-line entry point for the per-cell GWAS enrichment pipeline. This
#' is a thin wrapper around RunCellGWASEnrichment() (see R/run_enrichment_pipeline.R
#' for full parameter documentation) -- use this script directly from the
#' shell, cron, or a Docker CMD instead of editing example_run.R by hand.
#'
#' Usage:
#'   Rscript run_pipeline_cli.R --help
#'   Rscript run_pipeline_cli.R \
#'     --atac-rds data/my_atac_object.rds \
#'     --sumstats data/my_trait_sumstats.tsv.gz \
#'     --beta-col BETA --se-col SE --maf-col EAF --genome-build hg38 \
#'     --n-perm 200 --n-cores 4 \
#'     --output results/per_cell_gwas_enrichment.tsv

suppressPackageStartupMessages(library(optparse))

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
}, error = function(e) getwd())

source(file.path(script_dir, "R/utils_gwas.R"))
source(file.path(script_dir, "R/utils_fragments.R"))
source(file.path(script_dir, "R/utils_permutation.R"))
source(file.path(script_dir, "R/enrichment_engine.R"))
source(file.path(script_dir, "R/run_enrichment_pipeline.R"))
source(file.path(script_dir, "R/summarize_results.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

## ---- CLI options -----------------------------------------------------------
option_list <- list(
  make_option("--atac-rds", type = "character", default = NULL,
              help = "Path to a Seurat/Signac .rds object with peaks + attached fragments [required]"),
  make_option("--sumstats", type = "character", default = NULL,
              help = "Path to GWAS summary statistics (tab-delimited, optionally gzipped) [required]"),
  make_option("--assay", type = "character", default = NULL,
              help = "Assay name in the Seurat object [default: DefaultAssay(object)]"),
  make_option("--fragment-path", type = "character", default = NULL,
              help = paste("Optional override for the fragment file path(s), comma-separated if",
                            "multiple. Use this if the object's stored Fragment path (e.g. an HPC/NFS",
                            "path from wherever it was created) doesn't resolve inside the container",
                            "-- point it at your local fragments.tsv.gz (+ .tbi) instead.")),

  # GWAS sumstats column mapping
  make_option("--chr-col", type = "character", default = "CHR", help = "[default %default]"),
  make_option("--pos-col", type = "character", default = "POS", help = "[default %default]"),
  make_option("--beta-col", type = "character", default = "BETA", help = "[default %default]"),
  make_option("--se-col", type = "character", default = "SE", help = "[default %default]"),
  make_option("--z-col", type = "character", default = NULL,
              help = "Z-statistic column, used if beta/se are unavailable"),
  make_option("--pval-col", type = "character", default = NULL,
              help = "p-value column, used if beta/se/z are unavailable"),
  make_option("--maf-col", type = "character", default = NULL,
              help = "Allele frequency column (recommended)"),
  make_option("--snp-col", type = "character", default = NULL, help = "rsID column"),
  make_option("--genome-build", type = "character", default = "hg38",
              help = "hg38 or hg19; must match the ATAC object [default %default]"),
  make_option("--variant-panel-bed", type = "character", default = NULL,
              help = "Optional BED file to restrict GWAS variants to (e.g. HapMap3)"),
  make_option("--use-maf-covariate", action = "store_true", default = FALSE,
              help = "Include MAF (from --maf-col) as a regression covariate"),

  # Universe / peak restriction
  make_option("--whole-genome", action = "store_true", default = FALSE,
              help = "Use a whole-genome background instead of restricting to peaks (requires --chrom-sizes)"),
  make_option("--peaks-bed", type = "character", default = NULL,
              help = "Optional BED file overriding the object's own peaks (restrict_to_peaks mode only)"),
  make_option("--chrom-sizes", type = "character", default = NULL,
              help = "Two-column (chr, length) file, required with --whole-genome"),
  make_option("--blacklist-bed", type = "character", default = NULL,
              help = "Optional BED of regions to exclude (--whole-genome mode only)"),
  make_option("--consensus-peaks-bed", type = "character", default = NULL,
              help = paste("Optional BED of a consensus/reference peak set, added as a covariate to",
                            "absorb the generic 'GWAS signal is elevated in accessible chromatin'",
                            "effect. Strongly recommended with --whole-genome.")),

  # Cell grouping
  make_option("--group-by", type = "character", default = NULL,
              help = "Optional metadata column name to pool cells by (pseudobulk mode, e.g. cell_type) instead of per-cell"),
  make_option("--cells-file", type = "character", default = NULL,
              help = "Optional file with one cell barcode per line, to test a subset of cells (per-cell mode only)"),

  # Test parameters
  make_option("--n-perm", type = "integer", default = 200, help = "[default %default]"),
  make_option("--perm-method", type = "character", default = "circular",
              help = "circular or regioneR [default %default]"),
  make_option("--min-fragments", type = "integer", default = 200, help = "[default %default]"),
  make_option("--min-variants", type = "integer", default = 5, help = "[default %default]"),
  make_option("--n-cores", type = "integer", default = 1, help = "[default %default]"),
  make_option("--seed", type = "integer", default = 1, help = "[default %default]"),

  # Output
  make_option("--output", type = "character", default = "results/per_cell_gwas_enrichment.tsv",
              help = "[default %default]"),
  make_option("--plot-dir", type = "character", default = NULL,
              help = "Optional directory to write QC/summary plots (PDF) into"),
  make_option("--plot-group-by", type = "character", default = NULL,
              help = "Optional metadata/results column for the enrichment-by-group plot (e.g. cell_type)"),
  make_option("--quiet", action = "store_true", default = FALSE, help = "Suppress progress messages")
)

opt <- parse_args(OptionParser(option_list = option_list,
                                description = "Per-cell GWAS variant-effect-magnitude enrichment in ATAC fragments."))

if (is.null(opt[["atac-rds"]]) || is.null(opt$sumstats)) {
  stop("--atac-rds and --sumstats are required. Run with --help for full usage.", call. = FALSE)
}
if (opt[["whole-genome"]] && is.null(opt[["chrom-sizes"]])) {
  stop("--chrom-sizes is required when --whole-genome is set.", call. = FALSE)
}

verbose <- !opt$quiet

## ---- Load inputs ------------------------------------------------------------
if (verbose) message("Loading ATAC object: ", opt[["atac-rds"]])
atac <- readRDS(opt[["atac-rds"]])
assay <- if (is.null(opt$assay)) Seurat::DefaultAssay(atac) else opt$assay
Seurat::DefaultAssay(atac) <- assay

if (verbose) message("Loading GWAS summary statistics: ", opt$sumstats)
variant_panel <- if (!is.null(opt[["variant-panel-bed"]])) read_bed_as_granges(opt[["variant-panel-bed"]]) else NULL

gwas <- load_gwas_sumstats(
  sumstats_path = opt$sumstats,
  chr_col = opt[["chr-col"]], pos_col = opt[["pos-col"]],
  beta_col = opt[["beta-col"]], se_col = opt[["se-col"]],
  z_col = opt[["z-col"]], pval_col = opt[["pval-col"]],
  maf_col = opt[["maf-col"]], snp_col = opt[["snp-col"]],
  variant_panel = variant_panel,
  genome_build = opt[["genome-build"]]
)

covariates <- NULL
if (opt[["use-maf-covariate"]]) {
  if (is.null(opt[["maf-col"]])) stop("--use-maf-covariate requires --maf-col to be set.", call. = FALSE)
  covariates <- as.data.frame(GenomicRanges::mcols(gwas))["maf"]
}

peaks <- if (!is.null(opt[["peaks-bed"]])) read_bed_as_granges(opt[["peaks-bed"]]) else NULL
blacklist <- if (!is.null(opt[["blacklist-bed"]])) read_bed_as_granges(opt[["blacklist-bed"]]) else NULL
consensus_peaks <- if (!is.null(opt[["consensus-peaks-bed"]])) read_bed_as_granges(opt[["consensus-peaks-bed"]]) else NULL

chrom_sizes <- NULL
if (!is.null(opt[["chrom-sizes"]])) {
  chrom_sizes <- data.table::fread(opt[["chrom-sizes"]], header = FALSE, col.names = c("chr", "length"))
}

if (!is.null(opt[["group-by"]]) && !is.null(opt[["cells-file"]])) {
  warning("Both --group-by and --cells-file were supplied; --cells-file will be ignored ",
          "since --group-by (pseudobulk mode) takes precedence.")
}

cells <- NULL
if (!is.null(opt[["group-by"]])) {
  meta_col <- opt[["group-by"]]
  if (!meta_col %in% colnames(atac@meta.data)) {
    stop(sprintf("--group-by column '%s' not found in object metadata.", meta_col), call. = FALSE)
  }
  cells <- split(colnames(atac), atac@meta.data[[meta_col]])
  if (verbose) message(sprintf("Pooling cells into %d groups by '%s' (pseudobulk mode)", length(cells), meta_col))
} else if (!is.null(opt[["cells-file"]])) {
  cells <- readLines(opt[["cells-file"]])
  cells <- cells[nzchar(cells)]
  if (verbose) message(sprintf("Testing %d cells from --cells-file", length(cells)))
}

fragment_path <- NULL
if (!is.null(opt[["fragment-path"]])) {
  fragment_path <- trimws(strsplit(opt[["fragment-path"]], ",")[[1]])
  if (verbose) message("Overriding fragment path(s): ", paste(fragment_path, collapse = ", "))
}

## ---- Run the pipeline --------------------------------------------------------
res <- RunCellGWASEnrichment(
  object = atac,
  gwas_gr = gwas,
  fragment_path = fragment_path,
  assay = assay,
  restrict_to_peaks = !opt[["whole-genome"]],
  peaks = peaks,
  cells = cells,
  n_perm = opt[["n-perm"]],
  perm_method = opt[["perm-method"]],
  covariates = covariates,
  consensus_peaks = consensus_peaks,
  chrom_sizes = chrom_sizes,
  blacklist = blacklist,
  min_fragments = opt[["min-fragments"]],
  min_variants_in_annotation = opt[["min-variants"]],
  n_cores = opt[["n-cores"]],
  seed = opt$seed,
  verbose = verbose
)

## ---- Summarize + write outputs -----------------------------------------------
SummarizeCellGWASEnrichment(res$results)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
write.table(res$results, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote results to ", opt$output)

if (!is.null(opt[["plot-dir"]])) {
  dir.create(opt[["plot-dir"]], recursive = TRUE, showWarnings = FALSE)
  plots <- PlotCellGWASEnrichment(res$results, group_by = opt[["plot-group-by"]])
  for (nm in names(plots)) {
    out_pdf <- file.path(opt[["plot-dir"]], paste0(nm, ".pdf"))
    ggplot2::ggsave(out_pdf, plots[[nm]], width = 7, height = 5)
    message("Wrote plot to ", out_pdf)
  }
}
