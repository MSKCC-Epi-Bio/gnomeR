#' oncoKB_annotate
#' OncoKB annotator for MAF files.
#' @param maf A MAF file to be annotated.
#' @param cancer_types Data frame with samples mapped to cancer type for accurate levels of actionability. Default is NULL.
#' @param parallelize Should this be done in parallel. Default is TRUE
#' @return OncoKb annotate MAF.
#' @export
#' @examples
#' library(gnomeR)
#' maf <- mut[1:10,]
#' onco_maf <- oncoKB_annotate(maf)
#' @import
#' annotateMaf

oncoKB_annotate <- function(maf, cancer_types, parallelize = T){

  maf = maf %>%
    dplyr::mutate(dplyr::across(where(is.factor),as.character))
  return(annotateMaf::oncokb_annotate_maf(maf, cancer_types = NULL, parallelize = TRUE))
}

