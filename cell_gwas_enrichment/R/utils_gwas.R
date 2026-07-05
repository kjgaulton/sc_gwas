#' utils_gwas.R
#'
#' Loading and preparing GWAS summary statistics for enrichment testing.
#' The "variant effect magnitude" statistic used throughout the pipeline is
#' chi-square (chi2 = (beta/se)^2 = z^2), matching the convention used by
#' stratified LD score regression and MAGMA. Larger chi2 = stronger trait
#' association, regardless of effect direction.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
})

#' Load and prepare GWAS summary statistics as a GRanges object
#'
#' @param sumstats_path Path to a (optionally gzipped) tab/whitespace-delimited
#'   GWAS summary statistics file with a header row.
#' @param chr_col,pos_col Column names giving chromosome and base-pair position.
#'   Chromosome values are coerced to "chr"-prefixed UCSC style to match Signac
#'   fragment files (e.g. "1" -> "chr1").
#' @param beta_col,se_col Column names for effect size and its standard error.
#'   Used to compute chi2 = (beta/se)^2. Preferred input if available.
#' @param z_col Column name for a Z statistic, used to compute chi2 = z^2 if
#'   beta/se are not available.
#' @param pval_col Column name for a two-sided p-value, used to back out
#'   chi2 = qchisq(p, df = 1, lower.tail = FALSE) if beta/se/z are unavailable.
#'   This loses sign information, which is fine since chi2 is signless anyway,
#'   but is less numerically stable for very small p-values.
#' @param maf_col Optional column name for minor/effect allele frequency. If
#'   supplied it is carried through as a GRanges metadata column so it can be
#'   used later as a regression covariate (recommended: MAF is a known
#'   confounder of GWAS effect-size magnitude).
#' @param snp_col Optional rsID/variant-id column, kept as metadata for
#'   traceability.
#' @param variant_panel Optional character vector of variant IDs (matching
#'   snp_col) or a GRanges to restrict to (e.g. HapMap3 SNPs). Strongly
#'   recommended for genome-wide sumstats: restricting to a curated,
#'   well-imputed common-variant panel is standard practice for LDSC/MAGMA-
#'   style analyses and keeps the per-cell regression computationally
#'   tractable (see README, "Choosing a variant panel").
#' @param genome_build "hg38" or "hg19". Must match the genome build of the
#'   single-cell ATAC fragments -- the pipeline does NOT perform liftover.
#'
#' @return A GRanges of 1-bp variant positions with metadata columns:
#'   \code{chi2} (numeric), and optionally \code{maf}, \code{snp}.
load_gwas_sumstats <- function(sumstats_path,
                                chr_col = "CHR",
                                pos_col = "POS",
                                beta_col = "BETA",
                                se_col = "SE",
                                z_col = NULL,
                                pval_col = NULL,
                                maf_col = NULL,
                                snp_col = NULL,
                                variant_panel = NULL,
                                genome_build = c("hg38", "hg19")) {

  genome_build <- match.arg(genome_build)
  dt <- data.table::fread(sumstats_path, showProgress = FALSE)

  stopifnot(chr_col %in% names(dt), pos_col %in% names(dt))

  chi2 <- NULL
  if (beta_col %in% names(dt) && se_col %in% names(dt)) {
    chi2 <- (dt[[beta_col]] / dt[[se_col]])^2
  } else if (!is.null(z_col) && z_col %in% names(dt)) {
    chi2 <- dt[[z_col]]^2
  } else if (!is.null(pval_col) && pval_col %in% names(dt)) {
    chi2 <- stats::qchisq(dt[[pval_col]], df = 1, lower.tail = FALSE)
  } else {
    stop("Could not compute variant effect magnitude: provide beta_col+se_col, ",
         "z_col, or pval_col that exist in the sumstats file.")
  }

  chrom <- as.character(dt[[chr_col]])
  chrom <- ifelse(grepl("^chr", chrom, ignore.case = TRUE), chrom, paste0("chr", chrom))

  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(start = as.integer(dt[[pos_col]]), width = 1L)
  )
  GenomicRanges::mcols(gr)$chi2 <- chi2
  if (!is.null(maf_col) && maf_col %in% names(dt)) {
    GenomicRanges::mcols(gr)$maf <- dt[[maf_col]]
  }
  if (!is.null(snp_col) && snp_col %in% names(dt)) {
    GenomicRanges::mcols(gr)$snp <- dt[[snp_col]]
  }
  GenomeInfoDb::genome(gr) <- genome_build

  keep <- !is.na(gr$chi2) & is.finite(gr$chi2)
  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    message(sprintf("load_gwas_sumstats: dropping %d variants with missing/non-finite chi2", n_dropped))
  }
  gr <- gr[keep]

  if (!is.null(variant_panel)) {
    if (methods::is(variant_panel, "GRanges")) {
      gr <- IRanges::subsetByOverlaps(gr, variant_panel)
    } else if (!is.null(snp_col) && "snp" %in% names(GenomicRanges::mcols(gr))) {
      gr <- gr[gr$snp %in% variant_panel]
    } else {
      warning("variant_panel supplied but snp_col/snp metadata not available; ",
              "panel restriction skipped. Pass variant_panel as a GRanges to ",
              "restrict by position instead.")
    }
  }

  dup <- duplicated(gr)
  if (any(dup)) {
    message(sprintf("load_gwas_sumstats: removing %d duplicate-position variants", sum(dup)))
    gr <- gr[!dup]
  }

  message(sprintf("load_gwas_sumstats: %d variants loaded (build %s)", length(gr), genome_build))
  gr
}
