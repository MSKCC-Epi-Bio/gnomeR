#' Rename columns from API results to work with gnomeR functions
#'
#' Will return a named vector of internal column names as values and original data set names
#' as names as an attribute (`attr(x, "names_dict")`)
#' @param df_to_check A data frame to check and recode names as needed
#' @return a renamed data frame
#' @export
#' @examples
#'
#' rename_columns(df_to_check = gnomeR::mutations)
#' x <- rename_columns(df_to_check = gnomeR::sv)
#' attr(x, "names_dict")
rename_columns <- function(df_to_check) {

  names_df_long <- gnomeR::names_df %>%
    select(contains("_column_name")) %>%
    tidyr::pivot_longer(-"internal_column_name")


  which_to_replace <- intersect(names(df_to_check), unique(names_df_long$value))

  # create a temporary dictionary as a named vector- this should have all relevant values, including those unchanged
  names_dict <- names_df_long %>%
    dplyr::filter(.data$value %in% which_to_replace) %>%
    select("internal_column_name",  "value") %>%
    dplyr::distinct() %>%
    tibble::deframe()


  if(length(names_dict) > 0) {

    # store details on what has been changed.
    message <- purrr::map2_chr(names(names_dict),
                               names_dict,
                               ~paste0(.y, " renamed ", .x))

    names(message) <- rep("!", times = length(message))


    # rename those variables only
    df_to_check <- df_to_check %>%
      dplyr::rename(!!names_dict)

    attr(df_to_check, "names_dict") <- names_dict
  }

    return(df_to_check)
}




#' Utility Function to Extract SNV
#'
#' @param x string
#' @param n number of characters from right
#'
#' @return string
#' @noRd
#' @examples
#' substrRight("Hello", 2)
#'
substrRight <- function(x, n) {
  x <- as.character(x)
  substr(x, nchar(x) - n + 1, nchar(x))
}

#' Internal function to recode numeric CNA alteration values to factor values
#'
#' @param alteration_vector a vector of CNA alterations coded with any of the
#' following levels: neutral, deletion, amplification, gain, loss, homozygous deletion,
#' hemizygous deletion, loh, gain, high level amplification, 0, -1, -1.5, -2, 1, 2.
#'
#' @return a recoded CNA data set with factor alteration values. See details for code dictionary
#'
#' @details
#'
#' CNA is coded to the following key based on key: values below
#' - "neutral":  "0", "neutral",
#' - "deletion": "homozygous deletion", "-2",
#' - "deletion": "loh", "-1.5",
#' - "deletion": "hemizygous deletion", "-1",
#' - "amplification": "gain", "1",
#' - "amplification": high level amplification", "2",
#' @export
#' @examples
#' recode_cna(gnomeR::cna$alteration[1:10])

recode_cna <- function(alteration_vector){

  # *TODO - should we auto-recode Unknown/Unk to NA? - need to think on it

  # General Checks -------------------------------------------------------------
  # *TODO (I may remove these if this is just an internal function TBD )
  # (as it's already done in higher level function) but may change

  alteration_vector = as.character(alteration_vector)

  # CNA Levels Checks -----------------------------------------------------------
  levels_in_data <- tolower(names(table(alteration_vector)))

  # source: https://docs.cbioportal.org/file-formats/#data-file-1
  # python annotator ref with codes: https://github.com/oncokb/oncokb-annotator/blob/47e4a158ee843ead75445982532eb149db7f3106/AnnotatorCore.py#L158
   allowed_cna_levels <- tibble::tribble(
               ~detailed_coding, ~numeric_coding,   ~final_coding,
                      "neutral",          "0",        "neutral",
                    "deep loss",          "-2",      "deletion",
                    "deep loss",          "-1.5",    "deletion",
          "hemizygous deletion",          "-1",       "loss",
                         "gain",          "1",        "gain",
     "high level amplification",          "2",   "amplification")



  all_allowed <- unlist(allowed_cna_levels)
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }

  # Recode CNA Levels  ----------------------------------------------------------

  # create a named vector for recoding to final coding
  recode_values <- c(allowed_cna_levels$detailed_coding, allowed_cna_levels$numeric_coding)
  names(recode_values) <- c(allowed_cna_levels$final_coding, allowed_cna_levels$final_coding)

  recoded_alterations <- suppressWarnings(
    forcats::fct_recode(alteration_vector, !!!recode_values)
    )


    return(recoded_alterations)
  }




#' Create binary data.frames depending on type of mutation data
#'
#' @param data a dataset of alterations
#' @param samples a vector of unique sample ids
#' @param type a character indicator for which type of alteration the dataset contains
#' @return a data.frame of alterations
#' @keywords internal


.process_binary <- function(data,
                            samples,
                            type = c("mut", "del", "amp", "fus")){


  names_glue = switch(type,
                      mut =  rlang::expr("{hugo_symbol}"),
                      del = rlang::expr("{hugo_symbol}.Del"),
                      amp = rlang::expr("{hugo_symbol}.Amp"),
                      fus = rlang::expr("{hugo_symbol}.fus"))


  data_out <- data %>%
    filter(.data$sample_id %in% samples)

  data_out <- switch(type,
                     del = filter(data_out, .data$alteration %in% c("deletion")),
                     amp = filter(data_out, .data$alteration %in% c("amplification")),
                     mut = data_out,
                     fus = data_out)


  data_out %>%
    group_by(.data$sample_id,.data$hugo_symbol) %>%
    filter(row_number()==1) %>%
    mutate(fl = 1) %>%
    tidyr::pivot_wider(id_cols = "sample_id", names_from = "hugo_symbol", values_from  = "fl",
                       values_fill = 0, names_glue = rlang::eval_tidy(names_glue) ) %>%
    ungroup()
}


#' Check a Data Frame for Required Columns
#'
#' @param data A data frame to check
#' @param required_cols A character specifying names of columns to check
#' @param data_name Optionally specify how the data set should be called in error message.
#' Default is NULL and will call it a generic name.
#' @return If data set doesn't have required columns it will return an error message.
#' If it does have required columns, nothing will be returned
#' @keywords internal

.check_required_cols <- function(data, required_cols, data_name = NULL) {

  data_name <- data_name %||% ""
  column_names <- colnames(data)
  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your {data_name} data: {.field {which_missing}}")
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




