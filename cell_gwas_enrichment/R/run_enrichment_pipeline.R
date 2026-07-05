#' run_enrichment_pipeline.R
#'
#' Main driver: per-cell GWAS variant-effect-magnitude enrichment in ATAC
#' fragments, relative to a circular-permutation background.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(Seurat)
  library(Signac)
})

#' Run per-cell GWAS enrichment on ATAC fragments
#'
#' For every cell (or every group in \code{group_by}, see Details), tests
#' whether GWAS variant effect magnitude is enriched among variants
#' overlapping that cell's ATAC fragments, relative to a circular-permutation
#' null (see utils_permutation.R) fit with a regression model (see
#' enrichment_engine.R).
#'
#' @param object A Seurat object with a ChromatinAssay carrying attached
#'   Fragment file(s) (Signac::Fragments(object[[assay]])).
#' @param fragment_path Optional override for the fragment file path(s),
#'   bypassing whatever path(s) are stored on the object's Fragment
#'   object(s). Useful when the object was created on a different machine
#'   (e.g. an HPC/NFS path baked in at Fragment-object creation time that no
#'   longer resolves, such as inside a Docker container) -- point this at
#'   wherever you've actually placed the matching fragments.tsv.gz (+ .tbi
#'   index) instead of having to edit/resave the .rds. Accepts a single path
#'   or a character vector if the assay has multiple Fragment objects (e.g.
#'   one per sample). If NULL (default), paths are read from the object as
#'   usual via Signac::Fragments().
#' @param gwas_gr GRanges from load_gwas_sumstats(), with mcols$chi2 and
#'   optionally mcols$maf etc.
#' @param assay Assay name (default DefaultAssay(object)).
#' @param restrict_to_peaks Logical, default TRUE. If TRUE, both the variant
#'   universe and the permutation universe are restricted to the assay's
#'   peak set (Signac::granges(object[[assay]])), i.e. the test asks "are
#'   variants in THIS cell's accessible fragments more strongly associated
#'   than variants in OTHER accessible peaks not used by this cell". If
#'   FALSE, the universe is the whole genome (see chrom_sizes/blacklist),
#'   i.e. "accessible vs. inaccessible genome" -- this requires many more
#'   permutations and variants to be well powered and is much slower.
#' @param peaks Optional custom peak GRanges overriding the assay's own
#'   peaks (e.g. a pseudobulk peak set called across all cells).
#' @param cells Optional character vector of cell barcodes to test (default:
#'   all cells in object). Can also be a *named list* of barcode vectors, in
#'   which case each list element is treated as one pooled "cell" (its
#'   fragments are combined) -- this is how to run the identical test at
#'   pseudobulk/cluster/cell-type resolution instead of truly per-cell, for
#'   users who want more power at the cost of single-cell granularity.
#' @param n_perm Number of circular permutations per cell (default 200).
#'   The smallest attainable one-sided p-value is 1/(n_perm+1); increase for
#'   more resolution at the cost of runtime.
#' @param perm_method "circular" (default, fast) or "regioneR" (slower,
#'   independent randomization; requires the regioneR package).
#' @param covariates Optional data.frame/matrix aligned 1:1 with the rows of
#'   the \code{gwas_gr} argument as passed in (i.e. BEFORE universe
#'   restriction -- the function subsets it internally with the same mask
#'   used to restrict gwas_gr, so just pass e.g.
#'   \code{as.data.frame(mcols(gwas))["maf"]} using the same \code{gwas}
#'   object you pass as \code{gwas_gr}). Columns might include MAF and a
#'   local SNP-density LD proxy, partialled out of the regression as
#'   confounders.
#' @param consensus_peaks Optional GRanges of a consensus/reference peak set
#'   (e.g. the union of peaks across all cell types, a pseudobulk peak call,
#'   or an external atlas such as ENCODE cCREs). If supplied, a binary
#'   "in_consensus_peak" indicator (does variant j fall in ANY consensus
#'   peak) is added to the regression as an automatic covariate, in addition
#'   to anything passed via \code{covariates}.
#'
#'   This matters because GWAS signal is, in general, elevated in accessible
#'   chromatin/regulatory regions regardless of cell type -- so without this
#'   covariate, part of \code{obs_beta} for any given cell can simply reflect
#'   "variants in peaks are more GWAS-associated than variants outside any
#'   peak" rather than anything specific to that cell. The confound is worst
#'   with \code{restrict_to_peaks = FALSE}: a cell's real fragments are
#'   (by construction of ATAC-seq) almost entirely inside accessible
#'   chromatin, while a circular permutation across the *whole genome* will
#'   often land outside any peak, so nearly every cell can look "enriched"
#'   purely from the generic peaks-vs-genome effect. Including
#'   \code{consensus_peaks} lets the regression absorb that generic effect
#'   into its own coefficient, so \code{obs_beta} reflects enrichment
#'   specific to this cell's fragments over and above other cells'/generic
#'   accessible chromatin. Strongly recommended whenever
#'   \code{restrict_to_peaks = FALSE}; a warning is issued if it is omitted
#'   in that case. Under \code{restrict_to_peaks = TRUE} it is a no-op if
#'   \code{consensus_peaks} is (near-)identical to the peak universe itself
#'   (every tested variant is already "in a peak" by construction, so the
#'   indicator has ~zero variance) -- it's still useful there if you pass a
#'   *different*, broader/external consensus peak set than the one defining
#'   the universe.
#' @param chrom_sizes Required if restrict_to_peaks = FALSE: a data.frame
#'   with columns chr, length giving chromosome sizes for the genome build.
#' @param blacklist Optional GRanges of regions to exclude from the
#'   whole-genome universe (e.g. ENCODE blacklist), only used when
#'   restrict_to_peaks = FALSE.
#' @param min_fragments Minimum number of fragments a cell must have
#'   (post universe-restriction) to be tested (default 200).
#' @param min_variants_in_annotation Minimum number of GWAS variants that
#'   must overlap a cell's fragments to be tested (default 5). Cells below
#'   either threshold are reported with NA statistics and a `reason`.
#' @param n_cores Number of parallel workers (uses parallel::mclapply; set to
#'   1 to disable, which also works on Windows).
#' @param seed Random seed for reproducibility of permutations.
#' @param verbose Print progress messages.
#'
#' @return A list with:
#'   \item{results}{data.frame, one row per cell/group, with fragment counts,
#'     the enrichment statistics from fit_enrichment_regression(), an FDR
#'     (BH) q-value across all tested cells, and any cell metadata columns
#'     from object@meta.data.}
#'   \item{params}{list of parameters used, for provenance/reproducibility.}
RunCellGWASEnrichment <- function(object,
                                   gwas_gr,
                                   fragment_path = NULL,
                                   assay = Seurat::DefaultAssay(object),
                                   restrict_to_peaks = TRUE,
                                   peaks = NULL,
                                   cells = NULL,
                                   n_perm = 200,
                                   perm_method = c("circular", "regioneR"),
                                   covariates = NULL,
                                   consensus_peaks = NULL,
                                   chrom_sizes = NULL,
                                   blacklist = NULL,
                                   min_fragments = 200,
                                   min_variants_in_annotation = 5,
                                   n_cores = 1,
                                   seed = 1,
                                   verbose = TRUE) {

  perm_method <- match.arg(perm_method)
  set.seed(seed)

  ## 1. Define the universe (peaks or whole genome) ------------------------
  if (restrict_to_peaks) {
    universe <- if (!is.null(peaks)) reduce(peaks) else get_peak_universe(object, assay)
  } else {
    if (is.null(chrom_sizes)) {
      stop("chrom_sizes is required when restrict_to_peaks = FALSE (data.frame with chr, length).")
    }
    universe <- GenomicRanges::GRanges(chrom_sizes$chr,
                                        IRanges::IRanges(start = 1, end = chrom_sizes$length))
    if (!is.null(blacklist)) {
      universe <- setdiff(universe, reduce(blacklist))
    }
  }
  if (verbose) message(sprintf("Universe: %d regions, %s bp total",
                                length(universe), format(sum(as.numeric(width(universe))), big.mark = ",")))

  ## 1b. Reconcile chromosome-naming style (e.g. "chr1" vs "1") between the
  ## GWAS variants and the universe -- a very common mismatch that otherwise
  ## silently zeroes out every overlap regardless of genome build. No-op if
  ## they're already compatible.
  gwas_gr <- harmonize_seqnames(gwas_gr, universe, gr_label = "GWAS variants", reference_label = "ATAC peak/genome universe")

  ## 2. Restrict GWAS variants to the universe ------------------------------
  ## `covariates`, if supplied, must be aligned 1:1 with the *original* gwas_gr
  ## passed in (e.g. as.data.frame(mcols(gwas))["maf"]) -- it is subset here
  ## with the same logical mask used to restrict gwas_gr, so the two stay in
  ## lockstep regardless of how much the universe shrinks the variant set.
  keep <- IRanges::overlapsAny(gwas_gr, universe)
  if (!is.null(covariates)) {
    covariates_mat <- as.matrix(covariates)
    stopifnot(nrow(covariates_mat) == length(keep))
    covariates_mat <- covariates_mat[keep, , drop = FALSE]
  } else {
    covariates_mat <- NULL
  }
  gwas_gr <- gwas_gr[keep]
  if (length(gwas_gr) == 0) stop("No GWAS variants overlap the universe -- check genome build/chromosome naming.")
  if (verbose) message(sprintf("GWAS variants in universe: %d", length(gwas_gr)))

  ## 2b. Consensus-peak covariate: absorbs the generic "GWAS signal is
  ## elevated in accessible chromatin regardless of cell type" effect, so
  ## the per-cell coefficient reflects enrichment specific to that cell
  ## rather than to being in a peak at all. See the consensus_peaks roxygen
  ## docs above for why this matters most when restrict_to_peaks = FALSE.
  if (!is.null(consensus_peaks)) {
    in_consensus_peak <- as.numeric(IRanges::overlapsAny(gwas_gr, consensus_peaks))
    peak_var <- stats::var(in_consensus_peak)
    if (is.na(peak_var) || peak_var < 1e-8) {
      warning("consensus_peaks covariate has ~zero variance in the current universe ",
              "(nearly all or none of the tested variants fall inside it), so it will ",
              "act as a no-op once partialled against the intercept. This is expected if ",
              "restrict_to_peaks = TRUE and consensus_peaks matches the peak universe itself; ",
              "pass a different/broader consensus peak set, or use restrict_to_peaks = FALSE, ",
              "for this covariate to do anything.")
    }
    covariates_mat <- cbind(covariates_mat, in_consensus_peak = in_consensus_peak)
    if (verbose) {
      message(sprintf("Consensus-peak covariate: %d / %d variants (%.1f%%) in a consensus peak",
                       sum(in_consensus_peak), length(in_consensus_peak),
                       100 * mean(in_consensus_peak)))
    }
  } else if (!restrict_to_peaks && verbose) {
    warning("restrict_to_peaks = FALSE and consensus_peaks was not supplied: per-cell ",
            "enrichment estimates may largely reflect the generic 'GWAS signal is elevated ",
            "in accessible chromatin' effect rather than anything specific to individual ",
            "cells, since real fragments are almost always in peaks while genome-wide ",
            "permutations often are not. Consider passing consensus_peaks (e.g. the union ",
            "of peaks across all cells) to control for this.")
  }

  ## 3. Resolve cell groups (default: every cell individually) --------------
  if (is.null(cells)) {
    cell_groups <- as.list(colnames(object))
    names(cell_groups) <- colnames(object)
  } else if (is.list(cells)) {
    cell_groups <- cells
  } else {
    cell_groups <- as.list(cells)
    names(cell_groups) <- cells
  }
  all_barcodes <- unique(unlist(cell_groups))

  ## 4. Load + universe-clip fragments for all needed cells at once ---------
  if (verbose) message("Reading fragments...")
  frags_by_cell <- load_cell_fragments(object, universe, assay = assay, cells = all_barcodes,
                                        fragment_path = fragment_path)

  ## 5. Build the permutation universe map once, reused for every cell -----
  umap <- if (perm_method == "circular") build_universe_map(universe) else NULL

  ## 6. Per-group worker function --------------------------------------------
  chi2 <- gwas_gr$chi2
  run_one_group <- function(barcodes) {
    frag_list <- frags_by_cell[names(frags_by_cell) %in% barcodes]
    if (length(frag_list) == 0) {
      frag_gr <- GenomicRanges::GRanges()
    } else {
      frag_gr <- unlist(frag_list, use.names = FALSE)
    }
    n_frag <- length(frag_gr)

    if (n_frag < min_fragments) {
      return(data.frame(n_fragments = n_frag, n_variants_in_annotation = NA_integer_,
                         obs_beta = NA_real_, perm_mean = NA_real_, perm_sd = NA_real_,
                         z = NA_real_, p_empirical = NA_real_, p_empirical_two_sided = NA_real_,
                         enrichment_fold = NA_real_, n_perm = 0L,
                         reason = sprintf("fewer than min_fragments (%d < %d)", n_frag, min_fragments)))
    }

    annotation_obs <- IRanges::overlapsAny(gwas_gr, frag_gr)
    if (sum(annotation_obs) < min_variants_in_annotation) {
      return(data.frame(n_fragments = n_frag, n_variants_in_annotation = sum(annotation_obs),
                         obs_beta = NA_real_, perm_mean = NA_real_, perm_sd = NA_real_,
                         z = NA_real_, p_empirical = NA_real_, p_empirical_two_sided = NA_real_,
                         enrichment_fold = NA_real_, n_perm = 0L,
                         reason = sprintf("fewer than min_variants_in_annotation (%d < %d)",
                                           sum(annotation_obs), min_variants_in_annotation)))
    }

    perm_list <- if (perm_method == "circular") {
      generate_permutations(frag_gr, umap, n_perm)
    } else {
      generate_permutations_regioneR(frag_gr, universe, n_perm)
    }
    annotation_perm <- vapply(perm_list, function(pgr) IRanges::overlapsAny(gwas_gr, pgr),
                               FUN.VALUE = logical(length(gwas_gr)))
    colnames(annotation_perm) <- paste0("perm", seq_len(ncol(annotation_perm)))

    stats_row <- fit_enrichment_regression(chi2, annotation_obs, annotation_perm, covariates_mat)
    cbind(n_fragments = n_frag, stats_row, reason = NA_character_)
  }

  ## 7. Run over all groups (parallel optional) ------------------------------
  if (verbose) message(sprintf("Testing %d cell(s)/group(s)...", length(cell_groups)))
  if (n_cores > 1 && .Platform$OS.type == "unix") {
    rows <- parallel::mclapply(cell_groups, run_one_group, mc.cores = n_cores)
  } else {
    rows <- lapply(cell_groups, run_one_group)
  }

  results <- do.call(rbind, rows)
  results <- cbind(barcode = names(cell_groups), results, row.names = NULL)
  results$q_value <- stats::p.adjust(results$p_empirical, method = "BH")

  ## 8. Attach cell metadata, if this is truly per-cell (not pooled groups) -
  if (is.null(cells) || !is.list(cells)) {
    meta <- object@meta.data
    meta$barcode <- rownames(meta)
    results <- merge(results, meta, by = "barcode", all.x = TRUE, sort = FALSE)
  }

  list(
    results = results,
    params = list(assay = assay, restrict_to_peaks = restrict_to_peaks, n_perm = n_perm,
                  perm_method = perm_method, min_fragments = min_fragments,
                  min_variants_in_annotation = min_variants_in_annotation,
                  used_consensus_peak_covariate = !is.null(consensus_peaks),
                  n_variants_universe = length(gwas_gr), seed = seed,
                  date_run = Sys.time())
  )
}
