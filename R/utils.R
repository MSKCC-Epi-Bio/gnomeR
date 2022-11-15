#' Rename columns from API results to work with gnomeR functions
#'
#' @param df_to_check a data frame to check and recode names as needed
#'
#' @return a renamed data frame
#' @export
#' @examples
#'
#' rename_columns(df_to_check = gnomeR::mutations)
#' rename_columns(df_to_check = gnomeR::sv)
#'
rename_columns <- function(df_to_check) {

  names_df_long <- gnomeR::names_df %>%
    select(contains("_column_name")) %>%
    tidyr::pivot_longer(-"internal_column_name")


  which_to_replace <- intersect(names(df_to_check), unique(names_df_long$value))

  # create a temporary dictionary as a named vector
  temp_dict <- names_df_long %>%
    dplyr::filter(.data$value %in% which_to_replace) %>%
    select("internal_column_name",  "value") %>%
    dplyr::distinct() %>%
    tibble::deframe()


  if(length(temp_dict) > 0) {

    # store details on what has been changed.
    message <- purrr::map2_chr(names(temp_dict),
                               temp_dict,
                               ~paste0(.y, " renamed ", .x))

    names(message) <- rep("!", times = length(message))


    # rename those variables only
    df_to_check %>%
      dplyr::rename(!!temp_dict)
  }
}



#' Utility Function to Extract SNV
#'
#' @param x string
#' @param n number of characters from right
#'
#' @return string
#' @export
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
#' following levels: neutral, deletion, amplification, homozygous deletion,
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
  # python annotator ref: https://github.com/oncokb/oncokb-annotator/blob/47e4a158ee843ead75445982532eb149db7f3106/AnnotatorCore.py#L158
   allowed_cna_levels <- tibble::tribble(
               ~detailed_coding, ~numeric_coding,   ~final_coding,
                      "neutral",             "0",       "neutral",
          "homozygous deletion",            "-2",      "deletion",
                          "loh",          "-1.5",      "deletion",
          "hemizygous deletion",            "-1",      "deletion",
                         "gain",             "1", "amplification",
      "high level amplification",            "2", "amplification")


  # throw error if any values not in allowed list
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
    forcats::fct_recode(alteration_vector, !!!recode_values))

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
                     del = filter(data_out, .data$alteration %in% c("deletion","homozygous deletion","hemizygous deletion")),
                     amp = filter(data_out, .data$alteration %in% c("amplification","high level amplification")),
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



