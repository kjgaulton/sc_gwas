#' example_run.R
#'
#' End-to-end example. Replace the placeholder paths/columns with your own.

source("R/utils_gwas.R")
source("R/utils_fragments.R")
source("R/utils_permutation.R")
source("R/enrichment_engine.R")
source("R/run_enrichment_pipeline.R")
source("R/summarize_results.R")

library(Seurat)
library(Signac)

## --- 1. Load your Signac/Seurat ATAC object -------------------------------
## Must have a ChromatinAssay with peaks (granges(object)) and an attached
## Fragment object pointing at a tabix-indexed fragments.tsv.gz
## (Fragments(object[["peaks"]])).
atac <- readRDS("data/my_atac_object.rds")
DefaultAssay(atac) <- "peaks"

## --- 2. Load GWAS summary statistics ---------------------------------------
## Expects a header row; adjust column names to your file.
## Genome build MUST match the ATAC object (no liftover is performed here --
## use rtracklayer::liftOver first if they differ).
gwas <- load_gwas_sumstats(
  sumstats_path = "data/my_trait_sumstats.tsv.gz",
  chr_col = "CHR", pos_col = "POS",
  beta_col = "BETA", se_col = "SE",
  maf_col = "EAF", snp_col = "SNP",
  genome_build = "hg38"
)

## Recommended: restrict to a common, well-imputed variant panel (e.g.
## HapMap3) to keep the regression tractable and comparable to standard
## LDSC/MAGMA-style analyses.
## hm3 <- rtracklayer::import("data/hapmap3_variants.bed")
## gwas <- load_gwas_sumstats(..., variant_panel = hm3)

## --- 3. Run the enrichment test, per cell, restricted to peaks -------------
## Note on consensus_peaks (see RunCellGWASEnrichment() docs): it lets the
## regression absorb the generic "GWAS signal is elevated in accessible
## chromatin regardless of cell type" effect. It's NOT needed here, because
## restrict_to_peaks = TRUE already confines every tested variant to the
## peak universe -- passing the same peak set as consensus_peaks in this
## configuration would be a no-op (zero variance, since every variant is
## already "in a peak" by construction). It matters for the whole-genome
## variant (Section 3b below).
res <- RunCellGWASEnrichment(
  object = atac,
  gwas_gr = gwas,
  restrict_to_peaks = TRUE,      # set FALSE for a whole-genome background (slower, needs chrom_sizes)
  n_perm = 200,
  perm_method = "circular",
  covariates = as.data.frame(mcols(gwas))["maf"],  # optional; drop if maf not loaded
  min_fragments = 200,
  min_variants_in_annotation = 5,
  n_cores = 4,
  seed = 1
)

## --- 3b. Alternative: whole-genome background, controlling for peak status -
## Without restrict_to_peaks, real fragments are (by construction of
## ATAC-seq) almost always inside SOME peak, while genome-wide circular
## permutations often land outside any peak -- so nearly every cell would
## look "enriched" purely from the well-known generic peaks-vs-genome GWAS
## signal, not from anything specific to that cell. Passing consensus_peaks
## (here, the union of peaks across all cells) lets the regression absorb
## that generic effect so obs_beta isolates the cell-specific component.
##
## chrom_sizes <- data.frame(chr = paste0("chr", c(1:22, "X")),
##                            length = c(248956422, 242193529, ...))  # hg38
## res_genome <- RunCellGWASEnrichment(
##   object = atac, gwas_gr = gwas, restrict_to_peaks = FALSE,
##   chrom_sizes = chrom_sizes, consensus_peaks = get_peak_universe(atac),
##   n_perm = 1000, n_cores = 4
## )

## --- 4. Summarize + plot ----------------------------------------------------
SummarizeCellGWASEnrichment(res$results)

plots <- PlotCellGWASEnrichment(res$results, group_by = "cell_type")  # column in atac@meta.data
plots$p_hist
plots$p_group
plots$p_depth

write.table(res$results, "results/per_cell_gwas_enrichment.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

## --- Alternative: pseudobulk by cell type instead of per-cell -------------
## The same function supports pooling fragments across groups of cells by
## passing `cells` as a named list -- useful when per-cell power is too low
## and you want a more robust per-cell-type estimate instead.
##
## groups <- split(colnames(atac), atac$cell_type)
## res_pseudobulk <- RunCellGWASEnrichment(atac, gwas, cells = groups, n_perm = 1000)
