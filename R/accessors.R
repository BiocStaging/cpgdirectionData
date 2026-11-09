# Accessors for the cpgdirection data layers served from ExperimentHub.
#
# Two retrieval modes:
#   - ExperimentHub (default): the resource is fetched from the Hub on first
#     use and cached by ExperimentHub, so the download cost is paid once per
#     machine.
#   - local directory: options(cpgdirectionData.local = "<dir>") reads
#     "<dir>/<name>.rds" instead. For air-gapped machines, CI, and for
#     testing the full stack before the Hub upload exists.

# the resource registry: `name` is the stable identifier used by cpgdData()
# and by cpgdirection's data backend; `title` is the ExperimentHub title
# (it embeds the name, so pre-EH-ID queries by name resolve).
.CPGD_DATA_RESOURCES <- data.frame(
  name = c(
    "lookup_blood_hg19",
    "lookup_nasal_epithelium_hg19",
    "lookup_solid_tissue_hg19",
    "smr_directions",
    "cpg_positions_hg19",
    "brain_directions",
    "saliva_bridge_scores",
    "epicv2_probe_gene_annotation",
    "distal_links",
    "onco_eqtm_consensus"),
  title = c(
    "cpgdirection lookup_blood_hg19",
    "cpgdirection lookup_nasal_epithelium_hg19",
    "cpgdirection lookup_solid_tissue_hg19",
    "cpgdirection smr_directions",
    "cpgdirection cpg_positions_hg19",
    "cpgdirection brain_directions",
    "cpgdirection saliva_bridge_scores",
    "cpgdirection epicv2_probe_gene_annotation",
    "cpgdirection distal_links",
    "cpgdirection onco_eqtm_consensus"),
  description = c(
    paste("Blood catalogue lookup: 1,970,053 CpG-gene pairs (668,711 CpGs,",
          "hg19) with predicted methylation-expression direction, calibrated",
          "probability, confidence and evidence tier; one Ensembl gene per",
          "name at a strand-aware hg19 TSS, models retrained on harmonized",
          "pools (rebuilt in 0.99.7)"),
    paste("Nasal epithelium catalogue lookup: 1,500,942 CpG-gene pairs",
          "(552,167 CpGs, hg19; rebuilt in 0.99.7)"),
    paste("Solid tissue catalogue lookup: 1,184,630 CpG-gene pairs",
          "(541,536 CpGs, hg19; rebuilt in 0.99.7)"),
    paste("Two-sample SMR causal directions (GoDMC mQTL x eQTLGen cis-eQTL,",
          "blood): 413,982 CpG-gene pairs with tier, p_SMR, instruments and",
          "HEIDI results"),
    "EPIC v2 probe positions lifted to hg19: 930,178 probes",
    paste("Brain SMR causal directions (Brain-mMeta x BrainMeta v2):",
          "220,383 CpG-gene pairs with tier and score"),
    paste("Genome-wide synthesised saliva-bridge scores (periphery-to-brain",
          "concordance ensemble percentiles)"),
    paste("Multi-source EPIC v2 probe-to-gene annotation: 1,006,212 CpG-gene",
          "pairs (793,233 CpGs) as the union of five annotation tracks",
          "(Illumina UCSC RefGene, Illumina GencodeV41, Zhou GENCODE v41,",
          "Exeter GENCODE 49 (re-annotated manifest v3.0), Exeter regulatory",
          "elements) with per-track support flags, agreement counts and the",
          "Exeter Promoter_2000bp/Enhancer_5000bp labels"),
    paste("Distal CpG-to-gene links: ENCODE-rE2G thresholded enhancer-gene",
          "predictions (Gschwind et al. 2026) for 26 resting blood cell types",
          "and brain regions/cell types (47 released ENCODE annotation sets)",
          "plus scE2G single-cell enhancer-gene predictions (Sheth et al.",
          "2026) for nine PBMC cell types at E2G.Score.qnorm >= 0.177, all",
          "intersected with EPIC v2 probe positions (hg38) and keyed on the",
          "cg identifier: 10,205,122 rows, 1,289,742 CpG-gene pairs (525,222",
          "CpGs); one row per CpG x gene x source x biosample with the",
          "source score and a tissue class (blood/brain). A link proposes a",
          "target and carries no direction"),
    paste("Onco-eQTM cross-tumour cis-eQTM consensus (Korra et al. 2026;",
          "27 TCGA cancer types, 450K, hg19): 2,464,990 CpG-gene pairs with",
          "consensus direction, agreement tier (O1 >= 5 cancers agreeing:",
          "113,560; O2 2-4 agreeing; O3 disagreeing; O4 single cancer),",
          "number of cancer types and correlation summaries recovered from",
          "the t statistics. Reported beside, never inside, the peripheral",
          "direction call")),
  genome = c("hg19", "hg19", "hg19", "hg19", "hg19", "hg19", "hg19", "hg38",
             "hg38", "hg19"),
  stringsAsFactors = FALSE)


#' The data layers this package serves
#'
#' @return A \code{data.frame} with one row per resource: \code{name} (the
#'   identifier accepted by \code{\link{cpgdData}}), \code{title} (the
#'   ExperimentHub title), \code{description} and \code{genome}.
#' @examples
#' cpgdData_resources()
#' @export
cpgdData_resources <- function() {
  .CPGD_DATA_RESOURCES
}


#' Retrieve one cpgdirection data layer
#'
#' Fetches the named resource from ExperimentHub (cached locally by the Hub
#' after the first download), or from a local directory when
#' \code{options(cpgdirectionData.local = "<dir>")} is set -- the directory
#' must hold the resource as \code{"<name>.rds"}. The local mode serves
#' air-gapped machines and lets the cpgdirection stack run against a data
#' snapshot pinned for a specific analysis.
#'
#' @param name One of \code{cpgdData_resources()$name}.
#' @return A \code{data.frame}.
#' @examples
#' # The default retrieval is ExperimentHub (cached after first use):
#' #   pos <- cpgdData("cpg_positions_hg19")
#' # The local mode serves the same files from a directory, demonstrated
#' # here with a stand-in so the example runs without network access:
#' d <- tempfile(); dir.create(d)
#' saveRDS(data.frame(cpg_id = "cg00000029", chr = "chr16", pos = 53434200L),
#'         file.path(d, "cpg_positions_hg19.rds"))
#' old <- options(cpgdirectionData.local = d)
#' pos <- cpgdData("cpg_positions_hg19")
#' head(pos)
#' options(old)
#' @export
cpgdData <- function(name) {
  stopifnot(is.character(name), length(name) == 1L)
  reg <- .CPGD_DATA_RESOURCES
  if (!name %in% reg$name) {
    stop("Unknown resource '", name, "'. Available: ",
         paste(reg$name, collapse = ", "), call. = FALSE)
  }
  localdir <- getOption("cpgdirectionData.local",
                        Sys.getenv("CPGDIRECTIONDATA_LOCAL", ""))
  if (!is.null(localdir) && length(localdir) == 1L && nzchar(localdir)) {
    f <- file.path(localdir, paste0(name, ".rds"))
    if (!file.exists(f)) {
      stop("options(cpgdirectionData.local=) is set, but '", f,
           "' does not exist.", call. = FALSE)
    }
    return(readRDS(f))
  }
  eh <- ExperimentHub::ExperimentHub()
  hub <- AnnotationHub::query(eh, "cpgdirectionData")
  title <- reg$title[match(name, reg$name)]
  ids <- names(hub)[hub$title == title]
  if (!length(ids)) {
    stop("Resource '", name, "' is not listed on ExperimentHub for this ",
         "Bioconductor version (the Hub was reached and ", length(hub),
         " cpgdirectionData resource(s) are listed). For an offline copy, ",
         "point options(cpgdirectionData.local=) at a directory holding ",
         "the .rds files.", call. = FALSE)
  }
  # several versions of a resource can coexist on the Hub; take the newest
  eh[[ids[length(ids)]]]
}


#' @rdname layer-accessors
#' @name layer-accessors
#' @title Named accessors for the individual data layers
#' @description Convenience wrappers around \code{\link{cpgdData}}, one per
#'   resource. \code{getBloodLookup()}, \code{getNasalEpitheliumLookup()} and
#'   \code{getSolidTissueLookup()} return the per-tissue catalogue tables;
#'   \code{getSmrDirections()} the SMR causal directions;
#'   \code{getCpgPositionsHg19()} the EPIC v2 hg19 probe positions;
#'   \code{getBrainDirections()} the brain SMR directions;
#'   \code{getSalivaBridgeScores()} the saliva-bridge scores; and
#'   \code{getEpicv2ProbeGeneAnnotation()} the multi-source probe-to-gene
#'   annotation; \code{getDistalLinks()} the distal enhancer/loop
#'   CpG-to-gene links; and \code{getOncoEqtmConsensus()} the Onco-eQTM
#'   cross-tumour consensus.
#' @return A \code{data.frame}; see \code{\link{cpgdData_resources}} for the
#'   per-resource description.
#' @examples
#' # Hub retrieval (default):  pos <- getCpgPositionsHg19()
#' # Local-mode stand-in so the example runs without network access:
#' d <- tempfile(); dir.create(d)
#' saveRDS(data.frame(cpg_id = "cg00000029", chr = "chr16", pos = 53434200L),
#'         file.path(d, "cpg_positions_hg19.rds"))
#' old <- options(cpgdirectionData.local = d)
#' head(getCpgPositionsHg19())
#' options(old)
NULL

#' @rdname layer-accessors
#' @export
getBloodLookup <- function() cpgdData("lookup_blood_hg19")

#' @rdname layer-accessors
#' @export
getNasalEpitheliumLookup <- function() cpgdData("lookup_nasal_epithelium_hg19")

#' @rdname layer-accessors
#' @export
getSolidTissueLookup <- function() cpgdData("lookup_solid_tissue_hg19")

#' @rdname layer-accessors
#' @export
getSmrDirections <- function() cpgdData("smr_directions")

#' @rdname layer-accessors
#' @export
getCpgPositionsHg19 <- function() cpgdData("cpg_positions_hg19")

#' @rdname layer-accessors
#' @export
getBrainDirections <- function() cpgdData("brain_directions")

#' @rdname layer-accessors
#' @export
getSalivaBridgeScores <- function() cpgdData("saliva_bridge_scores")

#' @rdname layer-accessors
#' @export
getEpicv2ProbeGeneAnnotation <- function() cpgdData("epicv2_probe_gene_annotation")

#' @rdname layer-accessors
#' @export
getDistalLinks <- function() cpgdData("distal_links")

#' @rdname layer-accessors
#' @export
getOncoEqtmConsensus <- function() cpgdData("onco_eqtm_consensus")
