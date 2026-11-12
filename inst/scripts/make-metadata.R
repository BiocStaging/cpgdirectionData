# Regenerates inst/extdata/metadata.csv from the resource registry in
# R/accessors.R, so the two can never drift apart. Run from the package root:
#   Rscript inst/scripts/make-metadata.R <zenodo record id> <deposit version> \
#     [<name> ...]
# Without names, every row is written for that record (the Zenodo version
# holding all the .rds files). With names, only those resources are written
# for that record and every other row of the existing metadata.csv is kept
# as it is: a resource that did not change keeps the record it was
# registered from, so ExperimentHub is asked for the new resources only.
args <- commandArgs(trailingOnly = TRUE)
record <- if (length(args) >= 1) args[[1]] else "23257261"
deposit <- if (length(args) >= 2) args[[2]] else "2.99.8"
only <- if (length(args) >= 3) args[-(1:2)] else character(0)
stopifnot(grepl("^[0-9]+$", record), nzchar(deposit))
source("R/accessors.R")
reg <- .CPGD_DATA_RESOURCES
stopifnot(all(only %in% reg$name))
meta <- data.frame(
  Title = reg$title,
  Description = reg$description,
  BiocVersion = "3.23",
  Genome = reg$genome,
  SourceType = "CSV",
  SourceUrl = "https://doi.org/10.5281/zenodo.22225163",
  SourceVersion = deposit,
  Species = "Homo sapiens",
  TaxonomyId = 9606L,
  Coordinate_1_based = TRUE,
  DataProvider = "Ein-Dor lab, Reichman University (cpgdirection project)",
  Maintainer = "Tsachi Ein-Dor <teindor@runi.ac.il>",
  RDataClass = "data.frame",
  DispatchClass = "Rds",
  Location_Prefix = "https://zenodo.org/",
  RDataPath = paste0("records/", record, "/files/", reg$name, ".rds"),
  stringsAsFactors = FALSE)
if (length(only)) {
  f <- "inst/extdata/metadata.csv"
  stopifnot(file.exists(f))
  old <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = c(BiocVersion = "character",
                                 SourceVersion = "character"))
  stopifnot(identical(names(old), names(meta)))
  keep <- setdiff(reg$name, only)
  keep_title <- reg$title[match(keep, reg$name)]
  kept <- old[match(keep_title, old$Title), , drop = FALSE]
  if (anyNA(kept$Title)) {
    stop("metadata.csv has no row for: ",
         paste(keep[is.na(kept$Title)], collapse = ", "),
         ". Name it too, or write every row (no names).", call. = FALSE)
  }
  meta <- rbind(kept, meta[reg$name %in% only, , drop = FALSE])
  meta <- meta[match(reg$title, meta$Title), , drop = FALSE]
}
write.csv(meta, "inst/extdata/metadata.csv", row.names = FALSE)
what <- if (length(only)) paste(only, collapse = ", ") else "all"
cat("wrote inst/extdata/metadata.csv (", nrow(meta), " resources; ", what,
    " at record ", record, ", deposit ", deposit, ")\n", sep = "")
