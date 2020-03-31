#' oncoKB.annotate
#'
#' OncoKB annotator for MAF files.
#'
#' @param maf **NEED DESCRIPTION**
#' @return OncoKb annotate MAF.
#' @export
#'
#' @examples
#'
#' # NEED EXAMPLE, doesn't necessarily have to run

oncoKB.annotate <- function(input,output,lib.path = NULL){

  if(is.null(lib.path)) lib.path <- .libPaths()

  system(paste0("python ",lib.path,"/gnomeR/MafAnnotator.py -i ",input," -o ",output ))

}

