
# Recode Gene Aliases---------------------

#' Recode Hugo Symbol Column
#'
#'
#' Searches the Hugo Symbol column in a genomic dataframe to look for
#' any genes that have common gene name aliases,
#' and replaces those aliases with the accepted (most recent) gene name.
#' Function uses `gnomeR::impact_alias_table` by default as reference for
#' which aliases to replace and supports IMPACT panel alias replacement only at this time.
#' Custom tables can be provided as long as `hugo_symbol` and `alias` columns exist.
#'
#' @param genomic_df a gene_binary object
#' @param ... Other things passed
#'
#' @return A dataframe with a recoded Hugo Symbol columns
#' @export
#'
#' @examples
#' mut <- rename_columns(gnomeR::mutations)
#'
#' colnames(mut)
#'
#' colnames(recode_alias(genomic_df = mut))
#'
recode_alias <- function(genomic_df, alias_table = gnomeR::impact_alias_table, ...) {

  .check_required_cols(alias_table, "hugo_symbol", "alias")

  # get table of gene aliases (internal data)
  alias_table <- alias_table %>%
    dplyr::select("hugo_symbol", "alias")

  # recode aliases
  genomic_df$hugo_symbol_old <- genomic_df$hugo_symbol
  genomic_df$hugo_symbol <- purrr::map_chr(genomic_df$hugo_symbol,
                                           ~resolve_alias(gene_to_check = .x,
                                                          alias_table = alias_table))

  message <- genomic_df %>%
    dplyr::filter(.data$hugo_symbol_old != .data$hugo_symbol) %>%
    dplyr::select("hugo_symbol_old", "hugo_symbol") %>%
    dplyr::distinct()


  if(nrow(message) > 0) {
    vec_recode <- purrr::map2_chr(message$hugo_symbol_old,
                                  message$hugo_symbol,
                                  ~paste0(.x, " recoded to ", .y))

    names(vec_recode) <- rep("!", times = length(vec_recode))

    cli::cli_warn(c(
      "To ensure gene with multiple names/aliases are correctly grouped together, the
      following genes in your dataframe have been recoded (you can prevent this with {.code recode_aliases = FALSE}):",
      vec_recode))

  }

  genomic_df <- genomic_df %>%
    select(-"hugo_symbol_old")

  return(genomic_df)
}



#' Resolve Hugo Symbol Names with Aliases
#'
#' @param gene_to_check hugo_symbol to be check
#' @param alias_table table containing all the aliases
#'
#' @return if the accepted hugo symbol is input, it is returned back.
#' If an alias name is provided, the more common name/more up to date name is returned
#' @export
#'
#' @examples
#' resolve_alias("MLL4", alias_table = impact_alias_table)
#'
resolve_alias <- function(gene_to_check, alias_table) {

  if(gene_to_check %in% alias_table$alias) {

    alias_table %>%
      filter(.data$alias == gene_to_check) %>%
      pull("hugo_symbol") %>%
      first() %>%
      as.character()

  } else {
    as.character(gene_to_check)
  }
}

