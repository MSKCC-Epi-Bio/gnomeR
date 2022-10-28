# Check which panels gene is on----------------------

#' provide a list of impact panels a provided gene is found within
#'
#' @param hugo_symbol a vector of hugo symbols
#'
#' @return a data frame with hugo symbols and the IMPACT panels on which they are
#' included
#'
#' @examples
#'
#' hugos <- unique(gnomeR::mutations$hugoGeneSymbol)[1:10]
#'
#' which_impact_panel(hugos)
#'
#' @export

which_impact_panel <- function(hugo_symbol) {

  # get table of gene aliases (internal data)
  alias_table <- gnomeR::impact_alias_table %>%
    dplyr::select("hugo_symbol", "alias")

  # recode all genes to most common alias
  hugo_symbol <- purrr::map_chr(hugo_symbol,
                                ~resolve_alias(gene_to_check = .x,
                                               alias_table = alias_table))

  gene_panels <- gnomeR::gene_panels %>%
      filter(.data$gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505")) %>%
      select("gene_panel", "genes_in_panel") %>%
      tidyr::unnest(cols = c("genes_in_panel"))


  impact_results <- gene_panels %>%
    filter(.data$genes_in_panel %in% hugo_symbol) %>%
    distinct() %>%
    mutate(fill = "yes") %>%
    tidyr::pivot_wider(
      names_from = "gene_panel",
      values_from = "fill",
      values_fill = "no")

  impact_results

}
