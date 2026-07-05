#' utils_fragments.R
#'
#' Reading per-cell ATAC fragments from a Signac/Seurat object and (optionally)
#' restricting them to a peak universe. Fragments are read directly from the
#' underlying tabix-indexed fragments.tsv.gz file(s) rather than through
#' Signac's internal Fragment class, so this stays robust to Signac API
#' changes and works efficiently for very large fragment files (by only
#' pulling reads that fall in the peak universe via tabix, when restricting
#' to peaks).

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(Rsamtools)
})

#' Read a BED file (chrom, start, end, ...; extra columns ignored) into a
#' GRanges. Used by run_pipeline_cli.R to load peaks/blacklist/consensus-peak
#' arguments from plain BED files without requiring rtracklayer.
#'
#' @param path Path to a (optionally gzipped) BED file, no header.
#' @return GRanges, coordinates converted from BED's 0-based half-open to
#'   GRanges' 1-based inclusive.
read_bed_as_granges <- function(path) {
  dt <- data.table::fread(path, header = FALSE)
  GenomicRanges::GRanges(
    seqnames = as.character(dt[[1]]),
    ranges   = IRanges::IRanges(start = dt[[2]] + 1L, end = dt[[3]])
  )
}

#' Get the peak universe (accessible-chromatin regions) from a Signac object
#'
#' @param object A Seurat object with a ChromatinAssay.
#' @param assay Assay name (default: DefaultAssay(object)).
#' @return A reduced (merged, non-overlapping) GRanges of peaks.
get_peak_universe <- function(object, assay = Seurat::DefaultAssay(object)) {
  peaks <- Signac::granges(object[[assay]])
  reduce(peaks)
}

#' Locate the fragments.tsv.gz path(s) attached to a Signac assay
#'
#' @param object A Seurat object.
#' @param assay Assay name.
#' @return Character vector of fragment file paths.
get_fragment_paths <- function(object, assay = Seurat::DefaultAssay(object)) {
  frag_objs <- Signac::Fragments(object[[assay]])
  if (length(frag_objs) == 0) {
    stop("No Fragment objects found on assay '", assay, "'. Attach fragments with ",
         "Signac::CreateFragmentObject()/Fragments(object) <- ... before running the pipeline, ",
         "or pass fragment_path explicitly to RunCellGWASEnrichment().")
  }
  vapply(frag_objs, Signac::GetFragmentData, slot = "path", FUN.VALUE = character(1))
}

#' Read a fragments.tsv.gz file, optionally restricted to a set of regions
#'
#' @param fragment_path Path to a tabix-indexed fragments.tsv.gz file.
#' @param regions Optional GRanges; if supplied, only fragment records
#'   overlapping these regions are read (via tabix), which is dramatically
#'   faster/lighter than reading the whole file for peak-restricted analyses.
#' @param cells Optional character vector of cell barcodes to keep.
#' @return A data.table with columns chr, start, end, barcode, count
#'   (0-based, half-open, matching the 10x fragments.tsv.gz spec).
read_fragments <- function(fragment_path, regions = NULL, cells = NULL) {
  col_names <- c("chr", "start", "end", "barcode", "count")

  if (!is.null(regions)) {
    regions <- reduce(regions)
    tbx <- Rsamtools::TabixFile(fragment_path)
    open(tbx)
    on.exit(try(close(tbx), silent = TRUE))
    res <- Rsamtools::scanTabix(tbx, param = regions)
    lines <- unlist(res, use.names = FALSE)
    if (length(lines) == 0) {
      return(data.table::data.table(chr = character(), start = integer(),
                                     end = integer(), barcode = character(),
                                     count = integer()))
    }
    dt <- data.table::fread(text = lines, header = FALSE, col.names = col_names)
  } else {
    dt <- data.table::fread(fragment_path, header = FALSE, col.names = col_names)
  }

  if (!is.null(cells)) {
    dt <- dt[dt$barcode %in% cells, ]
  }
  dt
}

#' Convert a fragments data.table to a clipped, universe-restricted GRanges
#' and split it by cell barcode.
#'
#' @param frag_dt data.table from read_fragments() (0-based, half-open coords).
#' @param universe GRanges defining the allowed space fragments must live in
#'   (peak universe if restrict_to_peaks = TRUE, else the whole-genome/
#'   blacklist-masked universe). Fragments are clipped (via pintersect) to
#'   this universe so every fragment used downstream, including in the
#'   circular permutation, is fully contained in it.
#' @return A GRangesList of fragments, one element per cell barcode.
fragments_to_grangeslist <- function(frag_dt, universe) {
  if (nrow(frag_dt) == 0) return(GenomicRanges::GRangesList())

  # fragments.tsv.gz is 0-based half-open; convert to 1-based inclusive
  gr <- GenomicRanges::GRanges(
    seqnames = frag_dt$chr,
    ranges   = IRanges::IRanges(start = frag_dt$start + 1L, end = frag_dt$end)
  )
  S4Vectors::mcols(gr)$barcode <- frag_dt$barcode

  universe <- reduce(universe)
  hits <- GenomicRanges::findOverlaps(gr, universe)
  if (length(hits) == 0) return(GenomicRanges::GRangesList())

  clipped <- GenomicRanges::pintersect(
    gr[S4Vectors::queryHits(hits)],
    universe[S4Vectors::subjectHits(hits)]
  )
  S4Vectors::mcols(clipped) <- S4Vectors::mcols(gr[S4Vectors::queryHits(hits)])
  clipped <- clipped[width(clipped) > 0]

  split(clipped, clipped$barcode)
}

#' End-to-end helper: load fragments for an object and split by cell,
#' restricted to a universe.
#'
#' @param object Seurat/Signac object.
#' @param universe GRanges universe fragments must be clipped to.
#' @param assay Assay name.
#' @param cells Optional barcode subset (default: all cells in object).
#' @param fragment_path Optional override of the fragment file path(s).
#' @return GRangesList keyed by barcode.
load_cell_fragments <- function(object, universe, assay = Seurat::DefaultAssay(object),
                                 cells = NULL, fragment_path = NULL) {
  if (is.null(cells)) cells <- colnames(object)
  if (is.null(fragment_path)) fragment_path <- get_fragment_paths(object, assay)

  all_dt <- data.table::rbindlist(lapply(fragment_path, read_fragments,
                                          regions = universe, cells = cells))
  fragments_to_grangeslist(all_dt, universe)
}
