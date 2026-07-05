# cell_gwas_enrichment

Per-cell GWAS variant-effect-magnitude enrichment in single-cell ATAC-seq fragments, for a Signac/Seurat object.

## What this tests

For each cell, is the GWAS trait association signal (chi-square = (beta/se)^2, i.e. Z^2) elevated among variants that fall in that cell's accessible ATAC fragments, compared to a background where those same fragments have been placed randomly? This is a fragment-level analogue of stratified LD score regression / MAGMA-style annotation enrichment, run independently for every cell rather than for one annotation genome-wide.

Model fit per cell:

```
chi2_j = b0 + b1 * annotation_j + covariates_j'gamma + e_j
```

`annotation_j` = 1 if variant *j* overlaps the cell's fragments, else 0. `b1` is the enrichment coefficient. Significance is empirical: the same regression is refit on the same cell's fragments after permutation, many times, and `b1_observed` is ranked against that null. This avoids relying on OLS asymptotics for what is a small, discrete, LD-correlated annotation.

Optionally (default: on), both the variant universe and the permutation background are restricted to called ATAC peaks, so the test becomes "is this cell's used fraction of accessible chromatin more GWAS-enriched than the rest of accessible chromatin" rather than "accessible vs. inaccessible genome."

## Method details

- **Background/null**: circular permutation. A cell's whole fragment set is shifted together by a random offset within a condensed coordinate system built from the universe (peak set, or masked genome), wrapping at the boundary. Fragment widths, count, and relative spacing are exactly preserved — only their joint position changes. This is the same principle used by GoShifter (Trynka et al. 2015) for GWAS-annotation enrichment, and it preserves local LD structure between nearby fragments better than placing each fragment independently. An independent-randomization alternative (`regioneR::randomizeRegions`) is also provided (`perm_method = "regioneR"`), slower but a useful cross-check.
- **Regression engine**: covariates (if any — e.g. MAF, a local LD/SNP-density proxy) are partialled out of chi2 and out of every annotation column (observed + all permutations) in a single QR-based residualization (`stats::lm.fit`), then each permutation's coefficient is a simple residual-covariance ratio. This lets hundreds of permutations per cell be fit as one matrix operation instead of hundreds of `lm()` calls, and never materializes an n_variants x n_variants matrix.
- **Statistical validation**: `tests/validate_enrichment_stat.py` is a Python re-implementation of the same closed-form regression, checked against simulated data for (1) exact agreement with naive multivariate OLS, (2) uniform empirical p-values under the global null, (3) high power to detect an injected true signal. R itself was not available in the environment this pipeline was written in, so this validates the statistical method, not the R code verbatim — run it yourself (`python3 tests/validate_enrichment_stat.py`) and treat the R scripts as a straightforward, line-by-line port that should be smoke-tested against a small real or synthetic Signac object before a production run (see "Before running on real data" below).

## Files

```
R/utils_gwas.R              load_gwas_sumstats() -> GRanges with mcols$chi2
R/utils_fragments.R         get_peak_universe(), read_fragments(), load_cell_fragments()
R/utils_permutation.R       build_universe_map(), permute_fragments_circular(), generate_permutations()
R/enrichment_engine.R       fit_enrichment_regression() -- the core statistic
R/run_enrichment_pipeline.R RunCellGWASEnrichment() -- main driver
R/summarize_results.R       SummarizeCellGWASEnrichment(), PlotCellGWASEnrichment()
example_run.R               End-to-end example, written as R code you edit directly
run_pipeline_cli.R          Command-line wrapper around the same pipeline (see "Command-line usage")
tests/validate_enrichment_stat.py   Statistical calibration/power checks (Python, no R needed)
```

## Requirements

R packages: `Seurat`, `Signac`, `GenomicRanges`, `IRanges`, `S4Vectors`, `GenomeInfoDb`, `Rsamtools`, `data.table`, `ggplot2`. Optional: `regioneR` (only if using `perm_method = "regioneR"`), `parallel` (base R, for `n_cores > 1`), `optparse` (only needed for `run_pipeline_cli.R`).

```r
install.packages(c("Seurat", "Signac", "data.table", "ggplot2", "optparse"))
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors", "GenomeInfoDb", "Rsamtools"))
# optional:
BiocManager::install("regioneR")
```

## Command-line usage

`run_pipeline_cli.R` wraps the same `RunCellGWASEnrichment()` call as `example_run.R`, but takes everything as command-line flags instead of requiring you to edit R code — useful for cron/scheduler jobs, running inside the Docker container, or sweeping parameters from a shell loop.

```bash
Rscript run_pipeline_cli.R --help
```

Typical per-cell, peak-restricted run:

```bash
Rscript run_pipeline_cli.R \
  --atac-rds data/my_atac_object.rds \
  --sumstats data/my_trait_sumstats.tsv.gz \
  --chr-col CHR --pos-col POS --beta-col BETA --se-col SE \
  --maf-col EAF --use-maf-covariate \
  --genome-build hg38 \
  --n-perm 200 --n-cores 4 \
  --output results/per_cell_gwas_enrichment.tsv \
  --plot-dir results/plots --plot-group-by cell_type
```

Whole-genome background with the peaks-confound covariate (see "Controlling for the generic 'peaks are GWAS-enriched' effect" below):

```bash
Rscript run_pipeline_cli.R \
  --atac-rds data/my_atac_object.rds \
  --sumstats data/my_trait_sumstats.tsv.gz \
  --whole-genome --chrom-sizes data/hg38.chrom.sizes \
  --consensus-peaks-bed data/consensus_peaks.bed \
  --n-perm 1000 --n-cores 4 \
  --output results/per_cell_gwas_enrichment_genomewide.tsv
```

Pseudobulk by cell type instead of per-cell (pools fragments within each group before testing):

```bash
Rscript run_pipeline_cli.R \
  --atac-rds data/my_atac_object.rds \
  --sumstats data/my_trait_sumstats.tsv.gz \
  --group-by cell_type \
  --n-perm 1000 \
  --output results/pseudobulk_by_celltype.tsv
```

Peaks, blacklist, consensus-peaks, and variant-panel arguments all take plain BED files (chrom/start/end, no header); `--chrom-sizes` takes a two-column (chr, length) file like a standard UCSC `.chrom.sizes`. Run `--help` for the full flag list, which mirrors every argument documented on `RunCellGWASEnrichment()` in `R/run_enrichment_pipeline.R`.

## Running in Docker

A `Dockerfile` is included, built on `bioconductor/bioconductor_docker`, which comes with R + BiocManager preconfigured so the Bioconductor packages (GenomicRanges, Rsamtools, etc.) install cleanly. It also installs the system libraries (`libgsl-dev`, `libglpk-dev`, `libgmp-dev`, `libhdf5-dev`, `libxml2-dev`, `cmake`) that Seurat/Signac's compiled dependencies (`qlcMatrix`, `igraph`, `hdf5r`, ...) need to build from source, then Seurat and Signac themselves, each in its own layer with an immediate sanity check. Seurat is by far the slowest step — expect 15-20+ minutes on first build.

**Note:** this container was authored and reviewed for correctness, but not build-tested — the sandbox this was written in had no Docker daemon and no network access to CRAN/Bioconductor/Docker Hub. It has since been revised based on a real build failure (Signac failing to compile due to missing system libraries); if you hit another package install failure, the per-package sanity checks will tell you exactly which one, and you can debug it directly without waiting through a full rebuild:

```bash
# Comment out the failing package's later RUN steps in the Dockerfile so the
# build stops right after the layer that's failing, then:
docker build -t cell-gwas-debug .
docker run --rm -it cell-gwas-debug R -e 'install.packages("<package>", repos = "https://cloud.r-project.org")'
```

Earlier layers are cached, so this re-runs almost instantly up to the failing step and shows the real compiler/linker error instead of just "failed to install."

Build the image:

```bash
docker build -t cell-gwas-enrichment .
```

Run the example script, mounting local `data/` (your ATAC object + GWAS sumstats) and `results/` directories into the container:

```bash
docker run --rm -it \
  -v "$(pwd)/data:/pipeline/data" \
  -v "$(pwd)/results:/pipeline/results" \
  cell-gwas-enrichment Rscript example_run.R
```

Or drop into an interactive R session with the pipeline already sourced:

```bash
docker run --rm -it \
  -v "$(pwd)/data:/pipeline/data" \
  -v "$(pwd)/results:/pipeline/results" \
  cell-gwas-enrichment R
```

Or with `docker compose` (uses the volumes/command already configured in `docker-compose.yml`):

```bash
docker compose up --build
# or, for an interactive shell instead of running example_run.R:
docker compose run --rm cell-gwas-enrichment R
```

Place your Seurat/Signac `.rds` object and GWAS sumstats file under `./data` on the host before running — inside the container they'll appear at `/pipeline/data/...`, matching the paths used in `example_run.R`.

## Inputs

1. **Seurat/Signac object** with a `ChromatinAssay` that has peaks (`Signac::granges(object)`) and an attached `Fragment` object pointing at a tabix-indexed `fragments.tsv.gz` (`Signac::Fragments(object[["peaks"]])`). Standard 10x/Signac ATAC preprocessing output.
2. **GWAS summary statistics**, tab-delimited with a header, containing chromosome, position, and either (beta, se), a Z statistic, or a p-value — plus optionally MAF, used as a covariate. Genome build must match the ATAC object; this pipeline does **not** perform liftover.

## Usage

```r
source("R/utils_gwas.R"); source("R/utils_fragments.R"); source("R/utils_permutation.R")
source("R/enrichment_engine.R"); source("R/run_enrichment_pipeline.R"); source("R/summarize_results.R")

gwas <- load_gwas_sumstats("sumstats.tsv.gz", chr_col="CHR", pos_col="POS",
                            beta_col="BETA", se_col="SE", maf_col="EAF",
                            genome_build = "hg38")

res <- RunCellGWASEnrichment(atac_object, gwas, restrict_to_peaks = TRUE,
                              n_perm = 200, n_cores = 4)

SummarizeCellGWASEnrichment(res$results)
PlotCellGWASEnrichment(res$results, group_by = "cell_type")
```

See `example_run.R` for a fuller walkthrough, including how to instead run the identical test pooled by cluster/cell type (pass `cells = split(colnames(atac), atac$cell_type)`) if per-cell power turns out too low for your data.

## Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `restrict_to_peaks` | `TRUE` | Optional per the original request; `FALSE` requires `chrom_sizes` and is much slower/less powered (whole-genome universe) |
| `n_perm` | 200 | Smallest attainable one-sided p is `1/(n_perm+1)`; raise for more resolution |
| `perm_method` | `"circular"` | `"regioneR"` for slower independent-placement permutation |
| `covariates` | `NULL` | e.g. MAF; recommended if available |
| `consensus_peaks` | `NULL` | GRanges; absorbs the generic "GWAS signal is elevated in accessible chromatin regardless of cell type" effect (see below). Strongly recommended when `restrict_to_peaks = FALSE`. |
| `min_fragments` | 200 | Cells below this are skipped (NA + reason) |
| `min_variants_in_annotation` | 5 | Cells whose fragments overlap fewer GWAS variants than this are skipped |
| `n_cores` | 1 | Uses `parallel::mclapply` |

## Choosing a variant panel

Genome-wide sumstats can have 5-10M+ variants. Restricting to a curated, well-imputed, common-variant panel (e.g. HapMap3, ~1.2M SNPs — the standard choice for LDSC/MAGMA) via `load_gwas_sumstats(..., variant_panel = <GRanges or SNP ids>)` keeps the per-cell regression matrix a manageable size and avoids diluting the signal with poorly-imputed/rare variants.

## Controlling for the generic "peaks are GWAS-enriched" effect

GWAS signal is, in general, elevated in accessible chromatin/regulatory regions regardless of cell type — this is a well-established property of complex-trait genetics, not something specific to any one cell. Left unaddressed, this can dominate the per-cell test: with `restrict_to_peaks = TRUE`, both the tested variants and the permutation background already live inside the same peak universe, so this generic effect mostly cancels out. But with `restrict_to_peaks = FALSE` it does not — real fragments are, by construction of ATAC-seq, almost always inside some accessible region, while a genome-wide circular permutation will often land outside any peak. That alone can make nearly every cell look "enriched," independent of anything cell-specific.

`consensus_peaks` addresses this: pass a GRanges of a consensus/reference peak set (e.g. `get_peak_universe(object)` for the union of peaks across all cells, a pseudobulk peak call, or an external atlas like ENCODE cCREs), and the pipeline adds a binary "in_consensus_peak" indicator as an automatic regression covariate alongside anything passed via `covariates`. This lets the model absorb the generic peaks-vs-genome effect into its own coefficient, so `obs_beta` reflects enrichment specific to a cell's own fragments rather than just "these variants happen to be in accessible chromatin somewhere." A warning is printed if `restrict_to_peaks = FALSE` and `consensus_peaks` is omitted.

Under `restrict_to_peaks = TRUE`, this covariate is a no-op if `consensus_peaks` is (near-)identical to the peak universe itself — every tested variant is already "in a peak" by construction, so the indicator has ~zero variance and the pipeline warns you about this. It's still useful there if you supply a genuinely different/broader consensus set than the one defining the universe (e.g. an external reference peak atlas vs. this dataset's own peak calls).

## Statistical caveats (read before interpreting results)

- **Per-cell power is limited.** A single cell typically has thousands, not millions, of fragments, so the number of GWAS variants any one cell's fragments overlap can be small (tens to low hundreds even after peak restriction). Results for individual cells are noisy; treat per-cell p-values as exploratory, and look at aggregate patterns (e.g. `PlotCellGWASEnrichment(..., group_by = "cell_type")`) rather than single-cell calls. If you need well-powered single estimates, pool fragments by cluster/cell type using the `cells =` list argument.
- **LD is only partially accounted for.** Circular permutation preserves local spacing between a cell's own fragments but does not use an external LD reference panel, so it does not fully correct for GWAS signal clustering due to LD (the way stratified LDSC's LD-score regression does). Treat this as a first-pass, fragment-level enrichment test, not a substitute for a full baseline-model LDSC analysis if publication-grade heritability partitioning is the goal.
- **Depth confound.** Always check the `p_depth` plot from `PlotCellGWASEnrichment()` — spurious depth/enrichment correlation is a common artifact in per-cell analyses and should be reported/adjusted for if present (e.g. add a fragment-count covariate, or restrict comparisons to cells within similar depth bins).
- **Multiple testing** across potentially thousands of cells: `q_value` (BH FDR) is provided in `results`, but with `n_perm = 200` the p-value floor (0.005) limits how small an FDR-corrected q-value can get for any single cell; raise `n_perm` if you need to survive strict multiple-testing correction at single-cell resolution.
- **Genome build**: verify the ATAC fragments and GWAS sumstats use the same build; there is no automatic check beyond chromosome-name matching.

## Troubleshooting

- **"No GWAS variants overlap the universe" / a `.merge_two_Seqinfo_objects` warning about no shared sequence levels.** This is almost always a chromosome-naming mismatch (`"chr1"` vs `"1"`), not an actual genome-build mismatch -- `load_gwas_sumstats()` always produces `"chr"`-prefixed names, but ATAC peaks/fragments can just as easily use plain Ensembl/NCBI-style names depending on the reference used to build the object. `RunCellGWASEnrichment()` calls `harmonize_seqnames()` automatically to fix exactly this (it's a no-op if naming already matches), so if you still hit this error after updating, it likely means the two actually are on different genome builds (e.g. hg19 vs hg38) -- there's no automatic liftover, so you'd need to run one through `rtracklayer::liftOver()` first.
- **`TabixFile: file(s) do not exist` when reading fragments.** The `Fragment` object(s) attached to your Seurat/Signac object store an absolute path from wherever they were created (e.g. an HPC/NFS path), which won't resolve on a different machine or inside a Docker container. Either fix it permanently with `Signac::UpdatePath()` and resave the `.rds`, or use the `fragment_path` argument to `RunCellGWASEnrichment()` (`--fragment-path` on the CLI, comma-separated if there are multiple Fragment objects on the assay) to point at wherever you've actually placed the matching `fragments.tsv.gz` + `.tbi` index -- no need to touch the object itself.
- **Docker: editing an R script but the container still runs the old version.** The `Dockerfile` `COPY`s `R/` and `run_pipeline_cli.R` into the image at build time, so local edits need `docker build` again to take effect. For faster iteration, bind-mount them instead of relying on what's baked into the image: add `-v "$(pwd)/R:/pipeline/R" -v "$(pwd)/run_pipeline_cli.R:/pipeline/run_pipeline_cli.R"` to your `docker run` command.

## Before running on real data

Since this was developed without an R environment available to execute it end-to-end, smoke-test on your data (or a small synthetic Signac object with a handful of cells and a few thousand simulated variants) before a full production run, and skim through each `R/*.R` file — they're short and heavily commented.
