# CNA Recode -----------------------------------------------------


#' Function to recode numeric CNA alteration values to factor values
#'
#' @param alteration_vector a vector of CNA alterations coded with any of the
#' following levels: neutral, deletion, amplification, gain, loss, homozygous deletion,
#' hemizygous deletion, loh, gain, high level amplification, 0, -1, -1.5, -2, 1, 2.
#' @return a recoded CNA data set with factor alteration values. See details for code dictionary
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


# Rename Columns ----------------------------------------------------------


#' Rename columns from API results to work with gnomeR functions
#'
#' Will return a named vector of internal column names as values and original data set names
#' as names as an attribute (`attr(x, "names_dict")`)
#' @param df_to_check A data frame to check and recode names as needed
#' @return a renamed data frame
#' @export
#' @examples
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


# Extract Patient ID ------------------------------------------------------


#' Extract IMPACT Patient ID From Sample ID
#'
#' @param sample_id A character vector of IMPACT Tumor sample IDs
#'
#' @return Returns a vector of patient IDs
#' @export
#' @examples
#' sample_id = c("P-0000071-T01-IM3", "P-0000072-T02-IM4", "P-0000073-T03-IM5")
#' extract_patient_id(sample_id)
#'
extract_patient_id <- function(sample_id) {

  # Checks ----------------------------------------------------------------
  wrong_format <- sample_id[!stringr::str_detect(sample_id, "^P-\\d{1,}-T.*")]

  if (length(wrong_format) > 0) {
    cli::cli_abort("Some {.code sample_id} values do not match the expected IMPACT sample format (e.g `P-0000XX-T01-IM3`)")
  }

  patient_id = stringr::str_replace(sample_id, "-T.*", "")
  return(patient_id)
}


# Add point mutations ------------------------------------------------

#' Annotates point mutations of interest
#'
#' @param df Raw maf dataframe containing mutation data
#' @param gene_name Hugo symbol of gene of interest
#' @param chr_num Number of the chromosome with gene of interest
#' @param start_pos String referencing the start of the gene point location.
#' @param new_name String providing a name for this point mutation
#' @param mutually_exclusive Boolean determining if the point mutation should be
#' mutually exclusive from any other mutation on the gene or not. The default is `TRUE`.
#' @return  a data frame with updated hugo symbols for point mutations
#' @export
#'

add_point_mut <- function(df, gene_name, chr_num, start_pos, new_name,
                          mutually_exclusive = T){

  df <- .clean_and_check_cols(df) %>%
    mutate(start_position = as.character(start_position))

  if (mutually_exclusive){
    pm_df <- df %>%
      mutate(hugo_symbol = case_when(
        chromosome == chr_num & grepl(start_pos, start_position) ~ new_name,
        TRUE ~ hugo_symbol
      ))
  } else {

    # select only the rows that have the point mutation of interest
    pm_only <- df %>%
      mutate(hugo_symbol = case_when(
        chromosome == chr_num & grepl(start_pos, start_position) ~ new_name,
        TRUE ~ hugo_symbol
      )) %>%
      filter(.data$hugo_symbol == new_name)

    # The point mutation will be counted as both the original
    # hugo_symbol and also the new point mutation, so combine rows
    pm_df <- pf %>%
      rbind(pm_only)

  }

  return (pm_df)

}
