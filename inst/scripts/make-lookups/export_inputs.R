# Inputs of build_lookups.py that live in R packages, written as csv.gz.
#
#   Rscript export_inputs.R <out_dir> [<cpgdirection data dir>]
#
# Writes:
#   lookup_<tissue>.csv.gz   the shipped lookups (for the regression check)
#   epic_gene_chr.csv.gz     EPIC v2 annotation: probes per gene and chromosome
#   epic_zhou.csv.gz         EPIC v2 annotation: Zhou GENCODE v41 distance to TSS per probe-gene
#   measured_pairs.csv.gz    cpgd_measured_eqtms() with hg19 CpG positions
#   hgnc_table.csv.gz        HGNChelper::hgnc.table (symbol, approved symbol, chromosome)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "inputs"
if (length(args) >= 2) options(cpgdirection.data_dir = path.expand(args[[2]]))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("cpgdirection", quietly = TRUE), requireNamespace("HGNChelper", quietly = TRUE))
d <- function(b) as.data.table(cpgdirection:::.cpgd_data(b))
for (t in c("blood", "nasal_epithelium", "solid_tissue"))
  fwrite(d(paste0("lookup_", t, "_hg19")), file.path(out, paste0("lookup_", t, ".csv.gz")))
A <- d("epicv2_probe_gene_annotation")
fwrite(A[!is.na(chr) & nzchar(chr), .N, by = .(g = toupper(target_gene), chr)], file.path(out, "epic_gene_chr.csv.gz"))
fwrite(A[!is.na(zhou_dist_tss), .(cpg_id, g = toupper(target_gene), chr, zhou_dist_tss)], file.path(out, "epic_zhou.csv.gz"))
P <- d("cpg_positions_hg19")
M <- as.data.table(cpgdirection::cpgd_measured_eqtms())
M[, `:=`(chr = P$chr[match(cpg_id, P$cpg_id)], pos = P$pos[match(cpg_id, P$cpg_id)])]
fwrite(M[!is.na(pos), .(cpg_id, g = toupper(target_gene), chr, pos, tss_dist, tissue, cohorts)],
       file.path(out, "measured_pairs.csv.gz"))
utils::data("hgnc.table", package = "HGNChelper", envir = environment())
fwrite(hgnc.table, file.path(out, "hgnc_table.csv.gz"))
cat(sprintf("inputs written to %s (HGNChelper %s, cpgdirection %s)\n", out,
            as.character(utils::packageVersion("HGNChelper")), as.character(utils::packageVersion("cpgdirection"))))
