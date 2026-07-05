#' utils_permutation.R
#'
#' Generating a permuted (background) null for a cell's fragment set.
#'
#' Two methods are provided:
#'
#' 1. "circular" (default): a fast, vectorization-friendly circular shift of
#'    all of a cell's fragments together within a condensed "universe track"
#'    (the peak set, or the whole genome if restrict_to_peaks = FALSE).
#'    Fragment widths, count, and *relative spacing* are exactly preserved,
#'    only their joint offset changes. This is the same logic used by
#'    GoShifter (Trynka et al. 2015) for testing GWAS-variant enrichment in
#'    genomic annotations, and it has the advantage over independent random
#'    placement of better preserving local LD structure between nearby
#'    fragments. Implemented with plain vector arithmetic so it is cheap
#'    enough to call hundreds of times per cell.
#'
#' 2. "regioneR": wraps regioneR::randomizeRegions(), which independently
#'    re-places each fragment within the universe (with rejection sampling
#'    to keep them inside the mask). More conservative/independent than the
#'    circular shift, at the cost of speed. Optional dependency.
#'
#' Both return fragments fully contained within the supplied universe.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

#' Build a condensed coordinate map of a universe (peak set or genome)
#'
#' @param universe_gr GRanges, the allowed space for fragments (already
#'   reduced/merged is fine either way; this re-reduces internally).
#' @return A list describing the condensed 0-based coordinate system used by
#'   permute_fragments_circular().
build_universe_map <- function(universe_gr) {
  universe_gr <- sort(reduce(universe_gr))
  widths <- width(universe_gr)
  cum_end <- cumsum(as.numeric(widths))
  cum_start <- cum_end - widths
  list(universe = universe_gr, cum_start = cum_start, total_length = sum(as.numeric(widths)))
}

#' Circularly shift a cell's fragments within a universe map
#'
#' @param frag_gr GRanges of a single cell's fragments. Must already be fully
#'   contained within umap$universe (use fragments_to_grangeslist(), which
#'   clips to the universe, upstream).
#' @param umap Output of build_universe_map().
#' @return GRanges of permuted fragments, same widths/count as frag_gr.
permute_fragments_circular <- function(frag_gr, umap) {
  n <- length(frag_gr)
  if (n == 0) return(frag_gr)

  hits <- GenomicRanges::findOverlaps(frag_gr, umap$universe, type = "within", select = "first")
  if (any(is.na(hits))) {
    stop("permute_fragments_circular: some fragments are not fully contained in the ",
         "permutation universe. Clip fragments with fragments_to_grangeslist()/",
         "pintersect() against the same universe before permuting.")
  }

  u_start <- umap$cum_start[hits] + (start(frag_gr) - start(umap$universe)[hits])
  widths  <- width(frag_gr)

  offset <- stats::runif(1, min = 0, max = umap$total_length)
  shifted <- (u_start + offset) %% umap$total_length

  interval_idx <- findInterval(shifted, umap$cum_start, rightmost.closed = FALSE)
  interval_idx[interval_idx < 1] <- 1
  in_interval_offset <- shifted - umap$cum_start[interval_idx]

  new_start <- start(umap$universe)[interval_idx] + round(in_interval_offset)
  new_end   <- new_start + widths - 1L

  # Rare edge case: a shifted fragment overruns the end of the (typically
  # small) universe interval it landed in. Clip rather than let it spill
  # into an unrelated interval; this trims a small number of the widest
  # fragments landing in the narrowest peaks.
  interval_end <- end(umap$universe)[interval_idx]
  overflow <- pmax(new_end - interval_end, 0L)
  new_end <- new_end - overflow
  new_start <- pmin(new_start, new_end)

  GenomicRanges::GRanges(
    seqnames = seqnames(umap$universe)[interval_idx],
    ranges   = IRanges::IRanges(start = new_start, end = new_end)
  )
}

#' Generate n_perm circularly-permuted fragment sets for one cell
#'
#' @param frag_gr GRanges, one cell's (universe-clipped) fragments.
#' @param umap Output of build_universe_map(), built once and reused across
#'   all cells for efficiency.
#' @param n_perm Number of permutations.
#' @return A GRangesList of length n_perm.
generate_permutations <- function(frag_gr, umap, n_perm) {
  GenomicRanges::GRangesList(lapply(seq_len(n_perm), function(i) {
    permute_fragments_circular(frag_gr, umap)
  }))
}

#' Alternative: regioneR-based independent randomization (optional, slower)
#'
#' @param frag_gr GRanges, one cell's fragments.
#' @param universe_gr GRanges universe (peaks or genome minus blacklist).
#' @param n_perm Number of permutations.
#' @return A GRangesList of length n_perm.
generate_permutations_regioneR <- function(frag_gr, universe_gr, n_perm) {
  if (!requireNamespace("regioneR", quietly = TRUE)) {
    stop("regioneR is not installed. Install it (BiocManager::install('regioneR')) ",
         "or use perm_method = 'circular' instead.")
  }
  GenomicRanges::GRangesList(lapply(seq_len(n_perm), function(i) {
    regioneR::randomizeRegions(frag_gr, genome = universe_gr,
                                per.chromosome = TRUE, allow.overlaps = TRUE)
  }))
}
