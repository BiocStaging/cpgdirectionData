# Regenerates inst/extdata/metadata.csv from the resource registry in
# R/accessors.R, so the two can never drift apart. Run from the package root:
#   Rscript inst/scripts/make-metadata.R <zenodo record id> <deposit version>
# (both are those of the Zenodo version holding the ten .rds files).
args <- commandArgs(trailingOnly = TRUE)
record <- if (length(args) >= 1) args[[1]] else "23257261"
deposit <- if (length(args) >= 2) args[[2]] else "2.99.8"
stopifnot(grepl("^[0-9]+$", record), nzchar(deposit))
source("R/accessors.R")
reg <- .CPGD_DATA_RESOURCES
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
write.csv(meta, "inst/extdata/metadata.csv", row.names = FALSE)
cat("wrote inst/extdata/metadata.csv (", nrow(meta), " resources, record ", record, ", deposit ", deposit, ")\n", sep = "")
