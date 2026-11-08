# Write the rebuilt lookups as the cpgdirectionData Hub resources.
#   Rscript write_rds.R <build_lookups.py out dir> <hub_upload dir> [<shipped hub_upload dir>]
# Each lookup becomes <hub_upload>/lookup_<tissue>_hg19.rds: a data.frame with
# the shipped column names, order and classes, xz-compressed as make-data.R
# writes every resource. With a shipped directory, the column classes are
# checked against the shipped files.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
src <- args[[1]]; out <- args[[2]]; ref <- if (length(args) >= 3) args[[3]] else NULL
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cols <- c("cpg_id", "target_gene", "tss_dist", "status", "direction", "probability_plus1", "confidence", "evidence_tier")
for (t in c("blood", "nasal_epithelium", "solid_tissue")) {
  x <- fread(file.path(src, paste0("lookup_", t, ".csv.gz")), colClasses = list(character = c("cpg_id", "target_gene", "status", "evidence_tier")))
  x[, tss_dist := as.integer(tss_dist)]
  for (c in c("direction", "probability_plus1", "confidence")) set(x, j = c, value = as.numeric(x[[c]]))
  setcolorder(x, cols); x <- as.data.frame(x)
  if (!is.null(ref)) {
    s <- readRDS(file.path(ref, paste0("lookup_", t, "_hg19.rds")))
    stopifnot(identical(names(x), names(s)), identical(vapply(x, class, ""), vapply(s, class, "")))
  }
  f <- file.path(out, paste0("lookup_", t, "_hg19.rds"))
  saveRDS(x, f, compress = "xz")
  cat(sprintf("%-34s %9d rows  %6.1f MB  md5 %s\n", basename(f), nrow(x), file.size(f) / 1e6, unname(tools::md5sum(f))))
}
