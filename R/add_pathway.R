#' Pathway Alterations
#'
#' Checks if certain oncogenic signaling pathways are altered.  Pathways were curated from `https://pubmed.ncbi.nlm.nih.gov/29625050/`.  Please check for gene aliases before using.
#'
#' @param bin_mat a binary matrix
#' @param pathway The pathways to check.  The options are "RTK/RAS", "Nrf2", "PI3K", "TGFB", "p53", "Wnt", "Myc", "Cell cycle", "Hippo", "Notch", "all".
#'
#' @return a data frame: each sample is a row, columns are pathways, with values of 0/1 depending on pathway alteration status.
#' @export
#'
#' @examples
#'
#' binmat <- binary_matrix(mutation = mut, cna = cna, fusion = fusion)
#' pathway_df <- add_pathway(binmat, pathway = "all")
#'
add_pathway <- function(bin_mat, pathway = c("RTK/RAS", "Nrf2", "PI3K", "TGFB", "p53", "Wnt", "Myc", "Cell cycle", "Hippo", "Notch", "all")) {
  pathways <- gnomeR::pathways
  path_names <- names(pathways)

  if (any(str_detect(colnames(bin_mat), "\\.") == FALSE)) {
    colnames(bin_mat) <- c(paste0(colnames(bin_mat)[str_detect(colnames(bin_mat), "\\.", negate = TRUE)], ".mut"), colnames(bin_mat)[str_detect(colnames(bin_mat), "\\.")])
  }

  if (length(pathway) == 1) {
    if (pathway == "all") {
      path_out <- purrr::map_dfc(path_names, ~ .check_pathway(bin_mat, pathname = .x))
    }
  } else {
    if (length(pathway) > 1 & is.element("all", pathway)) {
      warning('"all" was included in the pathway list. "all" must be used alone. "all" was excluded, and pathway analysis was computed for the remaining pathways.')
      pathway <- pathway[str_detect(pathway, "all", negate = TRUE)]
      path_out <- purrr::map_dfc(pathway, ~ .check_pathway(bin_mat, pathname = .x))
    }

    path_out <- purrr::map_dfc(pathway, ~ .check_pathway(bin_mat, pathname = .x))
  }

  return(path_out)
}


.check_pathway <- function(bin_mat, pathname) {
  path_alt <- bin_mat %>%
    select(any_of(pathways[[pathname]])) %>%
    mutate(sum = rowSums(., na.rm = TRUE)) %>%
    transmute("{pathname}" := if_else(sum >= 1, 1, 0))

  return(path_alt)
}


