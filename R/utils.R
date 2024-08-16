


# Small Misc Utils  -----------------------------------------------------


#' Check a Data Frame for Required Columns
#'
#' @param data A data frame to check
#' @param required_cols A character specifying names of columns to check
#' @param add_to_message a vector (preferrably named) of text to add to the error message for specific cases
#' @return If data set doesn't have required columns it will return an error message.
#' If it does have required columns, nothing will be returned
#' @keywords internal

.check_required_cols <- function(data, required_cols, add_to_message = NULL) {

  # Get the name of the data object
  data_name <- deparse(substitute(data))

  column_names <- colnames(data)
  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    message <-
      c("Can't find required columns:", set_names(c(which_missing), "x"))

    add_to_message <- add_to_message %||% ""
    message <- c(message, add_to_message)
    cli::cli_abort(message)
  }
}



#' Add descriptive endings to hugo symbol names that do not have one already
#'
#' @param names hugo symbols to check
#' @param ending character ending to add to hugo symbol names without descriptive endings.
#' The default is ".mut". If interested in any type of alteration, use ".any".
#' @return a vector of hugo symbols where each entry has a descriptive ending
#' from the following list: ".Amp", ".Del", ".fus", ".cna", ".mut".
#' @keywords internal

.paste_endings = function(names, ending = NULL) {

  ending <- ending %||% ".mut"

  names[!str_detect(names, ".Amp|.Del|.fus|.cna")] <-
    paste0(stringr::str_trim(
      names[!str_detect(names, ".Amp|.Del|.fus|.cna")]), ending)

  return(names)
}

#' Check if all gene_binary columns except sample_id and other_vars are numeric
#'
#' @param alt_data a binary data frame created from `create_gene_binary()`
#' @return an error message if not all columns are numeric
#' @keywords internal

.abort_if_not_numeric = function(alt_data){

  # remove sample ID if it exists
  alt_data <- alt_data %>%
    select(-any_of("sample_id"))

  is_numeric <- purrr::map_lgl(alt_data, is.numeric)
  not_numeric <- names(is_numeric[!is_numeric])

  if(length(not_numeric) > 0) {
    cli::cli_abort("All alterations in your gene binary must be numeric and only can have values of 0, 1, or NA.
                   Please coerce the following columns to numeric or pass them to the `other_vars` argument before proceeding: {.field {not_numeric}}")
  }
}

