# The two inputs of build_pools.py that come from R packages, written as csv.gz.
#
#   Rscript export_inputs.R <out_dir> <annotables grch37.rda>
#
# Writes:
#   hgnc_table.csv.gz   HGNChelper::hgnc.table (symbol, approved symbol, chromosome)
#   ens_grch37.csv.gz   Ensembl GRCh37 genes from the annotables grch37 table (Ensembl Genes
#                       116, GRCh37): ensgene, symbol, chr, start, end, strand, strand-aware tss
#                       <https://raw.githubusercontent.com/stephenturner/annotables/master/data/grch37.rda>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: Rscript export_inputs.R <out_dir> <grch37.rda>" = length(args) >= 2 && file.exists(args[[2]]),
          requireNamespace("HGNChelper", quietly = TRUE))
out <- args[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::data("hgnc.table", package = "HGNChelper", envir = environment())
fwrite(hgnc.table, file.path(out, "hgnc_table.csv.gz"))
e <- new.env(); load(args[[2]], envir = e)
E <- as.data.table(e$grch37)[, .(ensgene, symbol, chr = paste0("chr", chr), start, end, strand,
                                 tss = ifelse(strand == 1, start, end))]
fwrite(unique(E), file.path(out, "ens_grch37.csv.gz"))
cat(sprintf("inputs written to %s (HGNChelper %s, annotables grch37 md5 %s)\n", out,
            as.character(utils::packageVersion("HGNChelper")), unname(tools::md5sum(args[[2]]))))
