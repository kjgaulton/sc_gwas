#' enrichment_engine.R
#'
#' Regression-based enrichment statistic (LDSC/MAGMA-style): for a given cell,
#' test whether GWAS variant effect magnitude (chi2) is elevated among
#' variants that overlap the cell's ATAC fragments ("in-annotation") relative
#' to variants that do not, using
#'
#'     chi2_j = b0 + b1 * annotation_j + covariates_j'gamma + e_j
#'
#' b1 is the enrichment coefficient of interest (excess chi2 per variant
#' attributable to being in an accessible fragment). Significance is assessed
#' empirically against a null distribution of b1 obtained by fitting the same
#' regression on many permuted (circularly-shifted) annotations for that same
#' cell, rather than relying on OLS asymptotics -- appropriate here because
#' annotations are small, discrete overlap indicators with unknown LD-induced
#' correlation structure among variants.
#'
#' All permutations for a cell are fit in a single matrix operation using the
#' Frisch-Waugh-Lovell (FWL) theorem: covariates are partialled out of chi2
#' and out of every annotation column (observed + permuted) via one QR-based
#' residualization (stats::lm.fit), and then each column's regression
#' coefficient reduces to a simple residual covariance ratio. This avoids
#' calling lm() hundreds of times per cell and avoids ever forming an
#' n_variants x n_variants matrix.

#' Fit the enrichment regression for one cell, observed + all permutations
#' at once.
#'
#' @param chi2 Numeric vector, length n_variants: GWAS variant effect
#'   magnitude for every variant in the (universe-restricted) analysis set.
#' @param annotation_obs Logical/0-1 numeric vector, length n_variants:
#'   whether each variant overlaps the cell's real fragments.
#' @param annotation_perm Logical/0-1 numeric matrix, n_variants x n_perm:
#'   whether each variant overlaps each permuted fragment set.
#' @param covariates Optional numeric matrix/data.frame, n_variants x n_cov
#'   (e.g. MAF, a local SNP-density LD proxy, distance to nearest TSS). Do not
#'   include an intercept column; one is added automatically.
#'
#' @return A one-row data.frame: obs_beta, perm_mean, perm_sd, z,
#'   p_empirical (one-sided, enrichment direction), p_empirical_two_sided,
#'   enrichment_fold, n_perm, n_variants_in_annotation.
fit_enrichment_regression <- function(chi2, annotation_obs, annotation_perm, covariates = NULL) {
  n <- length(chi2)
  stopifnot(length(annotation_obs) == n, nrow(annotation_perm) == n)

  annotation_obs <- as.numeric(annotation_obs)
  storage.mode(annotation_perm) <- "double"
  if (is.null(colnames(annotation_perm))) {
    colnames(annotation_perm) <- paste0("perm", seq_len(ncol(annotation_perm)))
  }
  X <- cbind(obs__ = annotation_obs, annotation_perm)

  if (is.null(covariates)) {
    C <- matrix(1, nrow = n, ncol = 1)
  } else {
    covariates <- as.matrix(covariates)
    C <- cbind(1, covariates)
  }

  fit <- stats::lm.fit(x = C, y = cbind(chi2__ = chi2, X))
  resid_all <- fit$residuals
  y_resid <- resid_all[, "chi2__"]
  X_resid <- resid_all[, colnames(X), drop = FALSE]

  num <- colSums(X_resid * y_resid)
  den <- colSums(X_resid^2)
  beta <- num / den

  obs_beta  <- unname(beta["obs__"])
  perm_beta <- beta[colnames(annotation_perm)]
  n_perm <- length(perm_beta)
  perm_beta <- perm_beta[is.finite(perm_beta)]
  if (length(perm_beta) < n_perm) {
    warning(sprintf("fit_enrichment_regression: dropped %d non-finite permutation betas",
                     n_perm - length(perm_beta)))
  }

  perm_mean <- mean(perm_beta)
  perm_sd   <- stats::sd(perm_beta)
  z <- (obs_beta - perm_mean) / perm_sd

  p_empirical           <- (1 + sum(perm_beta >= obs_beta)) / (length(perm_beta) + 1)
  p_empirical_two_sided <- (1 + sum(abs(perm_beta) >= abs(obs_beta))) / (length(perm_beta) + 1)

  n1 <- annotation_obs == 1
  enrichment_fold <- mean(chi2[n1]) / mean(chi2[!n1])

  data.frame(
    obs_beta = obs_beta,
    perm_mean = perm_mean,
    perm_sd = perm_sd,
    z = z,
    p_empirical = p_empirical,
    p_empirical_two_sided = p_empirical_two_sided,
    enrichment_fold = enrichment_fold,
    n_perm = length(perm_beta),
    n_variants_in_annotation = sum(n1)
  )
}
