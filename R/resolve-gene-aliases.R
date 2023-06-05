
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
#' @param alias_table a string indicating "impact", or a  dataframe with at least two columns (`hugo_symbol`,
#' `alias`) with one row for each pair.
#' @param supress_warnings If TRUE, function will return a list containing a dataframe of recoded results and a names vector of recoded aliases in data
#'
#' @return A dataframe with recoded Hugo Symbol columns
#' @export
#'
#' @examples
#' genomic_df <- rename_columns(gnomeR::mutations[1:5, ])
#'
#' alias_table <- data.frame("hugo_symbol" = c("New Symbol", "New Symbol2"),
#' "alias" = c("PARP1", "AKT1"))
#'
#' recode_alias(genomic_df, alias_table)
#' recode_alias(genomic_df, alias_table, supress_warnings = TRUE)

recode_alias <- function(genomic_df, alias_table = "impact", supress_warnings = FALSE) {

  # Checks ----------------------------------------------------

  .check_required_cols(genomic_df, "hugo_symbol")

  # make tibbles into data.frames
  if ("tbl" %in% class(alias_table)) {
    alias_table <- as.data.frame(alias_table)
  }

  alias_table <- switch(
    class(alias_table),
    "character" = {
      choices_arg <- c("impact", "IMPACT")
      lc = tolower(match.arg(alias_table, choices = choices_arg))
      switch(lc, "impact" = gnomeR::impact_alias_table)
      },

    "data.frame" = {
      .check_required_cols(alias_table, "hugo_symbol", "alias")
      alias_table
    })


  # make sure there is one gene per row
  if (is.character(alias_table$alias)) {
    if (any(stringr::str_detect(alias_table$alias, ","))) {
      cli::cli_abort(
        c("Error with {.code alias_table}. Are there multiple genes per row? You must provide a data frame with one gene-alias pair per row."),
        c("See {.code gnomeR::impact_alias_table} for an example on how to format data.")
      )
    }
  } else {
    cli::cli_abort("Error with {.code alias_table}. Did you provide a dataframe that has columns {.code hugo_symbol} and {.code alias}?")
  }

  # select only needed cols
  alias_table <- alias_table %>%
    dplyr::select("hugo_symbol", "alias")

  # Recode Aliases ---------------------
  genomic_df$hugo_symbol_old <- genomic_df$hugo_symbol
  genomic_df$hugo_symbol <- purrr::map_chr(genomic_df$hugo_symbol,
                                           ~resolve_alias(gene_to_check = .x,
                                                          alias_table = alias_table))

  changed_aliases <- genomic_df %>%
    dplyr::filter(.data$hugo_symbol_old != .data$hugo_symbol) %>%
    dplyr::select("hugo_symbol_old", "hugo_symbol") %>%
    dplyr::distinct()


  genomic_df <- genomic_df %>%
    select(-"hugo_symbol_old")

  result <- genomic_df

  # Aliases Warnings ------------------
  if (nrow(changed_aliases) > 0) {
    vec_recode <- purrr::map2_chr(
      changed_aliases$hugo_symbol_old,
      changed_aliases$hugo_symbol,
      ~ paste0(.x, " recoded to ", .y)
    )

    names(vec_recode) <- rep("!", times = length(vec_recode))
  } else {
    vec_recode <- c()
  }

  # format results according to supress warnings
  if (supress_warnings == TRUE) {
    result <- list("genomic_df" = genomic_df, "aliases_in_data" = vec_recode)
  }

  if (supress_warnings == FALSE) {
    if (nrow(changed_aliases) > 0) {
      cli::cli_warn(c(
        "To ensure gene with multiple names/aliases are correctly grouped together, the
        following genes in your dataframe have been recoded (if you are running {.code create_gene_binary()}
        you can prevent this with {.code alias_table = FALSE}):",
        vec_recode
      ))
    }
    result <- genomic_df
  }
  return(result)
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

