# Build the .rds resources served by cpgdirectionData.
#
# The three catalogue lookups (lookup_<tissue>_hg19) are not converted here
# since 0.99.5: inst/scripts/make-lookups/ rebuilds them from their sources
# and writes the .rds files directly (see its README.md). This script
# converts six of the other nine; distal_links and onco_eqtm_consensus are
# written by cpgdirection's tools/build_distal_links.R and
# tools/build_onco_eqtm.R, and smr_directions is frozen (see below).
#
# The layers are converted 1:1 from the csv.gz files shipped in the FULL
# (GitHub) build of the cpgdirection package -- themselves produced by the
# pipelines documented in the cpgdirection repository (catalogue models, SMR
# analyses, and the multi-source EPIC v2 annotation build that merges the
# Illumina EPIC-8v2-0_A2 manifest, the Zhou lab SeSAMe GENCODE v41
# annotation, and the Exeter GENCODE 47 re-annotation, Zenodo 15181885).
#
# Since 0.99.8 two more resources are converted here, both from cpgdirection
# 1.1.0 or later:
#   cpg_coordinates    <extdata_dir>/cpg_coordinates.csv.gz, built by
#                      tools/build_cpg_coordinates.R;
#   smr_directions_v2  the second SMR layer (three expression outcomes,
#                      replication rung, per-pair accuracy), built by
#                      tools/build_smr_layer_interval.py and named
#                      smr_directions.csv.gz on the cpgdirection side. It is
#                      not bundled with cpgdirection, so its path is the third
#                      argument.
# The first SMR table (smr_directions, 413,982 pairs) is frozen as deposited
# in Zenodo version 2.99.10: cpgdirection 0.99.x and 1.0.x take every row of
# that resource, so the second layer is served under its own name and
# smr_directions is not rebuilt.
#
# Usage:
#   Rscript inst/scripts/make-data.R <extdata_dir> <out_dir> [<smr_v2.csv.gz>]
# where <extdata_dir> is inst/extdata of a full cpgdirection checkout (or
# system.file("extdata", package = "cpgdirection") of a full install),
# <out_dir> receives the .rds files deposited on Zenodo, and <smr_v2.csv.gz>
# is the second SMR layer. A resource whose source file is absent is skipped
# and named.

args <- commandArgs(trailingOnly = TRUE)
extdata <- if (length(args) >= 1) args[[1]] else
  system.file("extdata", package = "cpgdirection")
out_dir <- if (length(args) >= 2) args[[2]] else "hub_upload"
smr_v2 <- if (length(args) >= 3) args[[3]] else ""
stopifnot(nzchar(extdata), dir.exists(extdata))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

resources <- c(
  cpg_positions_hg19 = file.path(extdata, "cpg_positions_hg19.csv.gz"),
  brain_directions = file.path(extdata, "brain_directions.csv.gz"),
  saliva_bridge_scores = file.path(extdata, "saliva_bridge_scores.csv.gz"),
  epicv2_probe_gene_annotation =
    file.path(extdata, "epicv2_probe_gene_annotation.csv.gz"),
  cpg_coordinates = file.path(extdata, "cpg_coordinates.csv.gz"),
  smr_directions_v2 = smr_v2)

for (r in names(resources)) {
  src <- resources[[r]]
  if (!nzchar(src) || !file.exists(src)) {
    cat(sprintf("%-32s skipped (no source file)\n", r))
    next
  }
  x <- as.data.frame(data.table::fread(src, showProgress = FALSE))
  out <- file.path(out_dir, paste0(r, ".rds"))
  saveRDS(x, out, compress = "xz")
  cat(sprintf("%-32s %8d x %2d  ->  %s (%.1f MB)\n",
              r, nrow(x), ncol(x), out, file.size(out) / 1e6))
}
cat("done. Deposit the contents of '", out_dir,
    "' as a new version of the Zenodo record, then run make-metadata.R.\n",
    sep = "")
