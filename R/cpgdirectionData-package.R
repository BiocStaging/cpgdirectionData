#' cpgdirectionData: data layers for the cpgdirection package
#'
#' ExperimentHub companion to \pkg{cpgdirection}. Bioconductor software
#' packages are capped at 5 MB per file; the twelve data layers served here
#' are larger than that, so they live on the Hub and are fetched on first
#' use (then cached locally by ExperimentHub). The \pkg{cpgdirection}
#' package retrieves them transparently through its data backend; direct
#' access goes through \code{\link{cpgdData}} or the named accessors.
#'
#' Provenance: the layers are built by the scripts under
#' \code{inst/scripts/} from the sources documented there (the cpgdirection
#' catalogues and SMR analyses; the Illumina EPIC-8v2-0_A2 manifest; the
#' Zhou lab SeSAMe EPIC v2 annotation and masks; the Exeter EPIC v2
#' re-annotation, Zenodo record 15181885; the Garvan EPICv2manifest
#' cross-hybridization evidence).
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @docType package
#' @name cpgdirectionData
#' @aliases cpgdirectionData-package
"_PACKAGE"
