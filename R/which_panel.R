# Check which panels gene is on----------------------

#' provide a list of impact panels a provided gene is found within
#'
#' @param genomic_df a data frame containing at least one column with hugo_symbols
#' @param impact_only indicator to only check IMPACT panels (default = T)
#'
#' @return a data frame with two columns: 1) the original list of hugo_symbols (recoded as more common aliases)
#' 2) a list of all the panels that hugo_symbol is found within
#'
#' @examples
#'
#' #select first 6 unique hugo symbols from example dataset
#' mut6 <- gnomeR::mutations %>%
#'   rename_columns()%>%
#'   select(hugo_symbol)%>%
#'   unique()%>%
#'   head()
#'
#' which_panel(mut6)
#' which_panel(mut6, impact_only = F)
#'
#' #example with some uncommon hugo_symbols
#' hugo_symbol <- c("ZNRF3", "MLL3")
#' example <- as.data.frame(hugo_symbol)
#'
#' which_panel(example)
#'
#' @export

which_panel <- function(genomic_df, impact_only = T) {
  # recode all genes to most common alias
  genomic_df <- recode_alias(genomic_df)%>%
    select("hugo_symbol")

  if (impact_only) {
    gene_panels <- gnomeR::gene_panels %>%
      filter(.data$gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505"))
  } else{
    gene_panels <- gnomeR::gene_panels
  }

  genomic_df <- genomic_df %>%
    mutate(panels = NA)

  # empty
  vec_panels <- c()

  for (gene in genomic_df$hugo_symbol) {
    for (panel_name in gene_panels$gene_panel) {
      if (gene %in% gene_panels$genes_in_panel[gene_panels$gene_panel == panel_name][[1]]) {
        vec_panels <- c(vec_panels, panel_name)
      }
    }

    if (length(vec_panels) > 0) {
      genomic_df$panels[genomic_df$hugo_symbol == gene] <- list(vec_panels)
    } else {
      genomic_df$panels[genomic_df$hugo_symbol == gene] <- NA
    }
  }

  return(genomic_df)
}
