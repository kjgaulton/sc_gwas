#' summarize_results.R
#'
#' Summary tables and QC/interpretation plots for RunCellGWASEnrichment()
#' output.

suppressPackageStartupMessages({
  library(ggplot2)
})

#' Print a text summary of per-cell enrichment results
#'
#' @param results The $results data.frame from RunCellGWASEnrichment().
#' @param alpha Significance threshold for the "nominally significant" count.
SummarizeCellGWASEnrichment <- function(results, alpha = 0.05) {
  tested <- results[!is.na(results$p_empirical), ]
  n_total <- nrow(results)
  n_tested <- nrow(tested)
  n_sig <- sum(tested$p_empirical < alpha)
  n_sig_fdr <- sum(tested$q_value < alpha, na.rm = TRUE)

  cat(sprintf("Cells/groups provided:        %d\n", n_total))
  cat(sprintf("Cells/groups tested:          %d (%.1f%%)\n", n_tested, 100 * n_tested / n_total))
  cat(sprintf("  -- skipped (fragments/variants below threshold): %d\n", n_total - n_tested))
  if (n_tested > 0) {
    cat(sprintf("Nominally enriched (p < %.2f):       %d (%.1f%% of tested)\n",
                alpha, n_sig, 100 * n_sig / n_tested))
    cat(sprintf("FDR-significant (q < %.2f):          %d (%.1f%% of tested)\n",
                alpha, n_sig_fdr, 100 * n_sig_fdr / n_tested))
    cat(sprintf("Median enrichment fold (tested cells): %.3f\n", stats::median(tested$enrichment_fold)))
    cat(sprintf("Median obs_beta (tested cells):        %.4g\n", stats::median(tested$obs_beta)))
  }
  invisible(tested)
}

#' QC/interpretation plots for per-cell enrichment results
#'
#' Produces:
#'  (1) a histogram of empirical p-values -- under no true enrichment signal
#'      anywhere, this should look roughly uniform; a pileup near 0 across
#'      many cells suggests a real, broadly shared enrichment signal (or a
#'      systematic artifact -- check plot 3 before concluding the former).
#'  (2) if group_by is supplied (e.g. cell type/cluster column in results),
#'      a boxplot of obs_beta by group.
#'  (3) a scatter of n_fragments vs. obs_beta, to check whether enrichment
#'      strength is confounded with sequencing depth per cell (a common
#'      artifact in per-cell analyses) -- there should be no strong trend.
#'
#' @param results The $results data.frame from RunCellGWASEnrichment().
#' @param group_by Optional column name in results to group plot (2) by.
#' @return A named list of ggplot objects: p_hist, p_group (if requested),
#'   p_depth.
PlotCellGWASEnrichment <- function(results, group_by = NULL) {
  tested <- results[!is.na(results$p_empirical), ]
  plots <- list()

  plots$p_hist <- ggplot2::ggplot(tested, ggplot2::aes(x = p_empirical)) +
    ggplot2::geom_histogram(bins = 30, boundary = 0, fill = "steelblue", color = "white") +
    ggplot2::labs(title = "Per-cell empirical enrichment p-values",
                  x = "p (empirical, one-sided)", y = "Number of cells") +
    ggplot2::theme_minimal()

  if (!is.null(group_by) && group_by %in% names(tested)) {
    plots$p_group <- ggplot2::ggplot(tested, ggplot2::aes(x = .data[[group_by]], y = obs_beta)) +
      ggplot2::geom_boxplot(outlier.size = 0.5) +
      ggplot2::labs(title = "Enrichment coefficient by group",
                    x = group_by, y = "obs_beta (regression enrichment coefficient)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  plots$p_depth <- ggplot2::ggplot(tested, ggplot2::aes(x = n_fragments, y = obs_beta)) +
    ggplot2::geom_point(alpha = 0.3, size = 0.7) +
    ggplot2::scale_x_log10() +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "firebrick") +
    ggplot2::labs(title = "Enrichment vs. sequencing depth (QC)",
                  x = "Fragments per cell (log10)", y = "obs_beta") +
    ggplot2::theme_minimal()

  plots
}
