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


#' Convert formula selector to a named list
#' Functions takes a list of formulas, a named list, or a combination of named
#' elements with formula elements and returns a named list.
#' For example, `list(age = 1, starts_with("stage") ~ 2)`.
#'
#' @section Shortcuts:
#' A shortcut for specifying an option be applied to all columns/variables
#' is omitting the LHS of the formula.
#' For example, `list(~ 1)` is equivalent to passing `list(everything() ~ 1)`.
#'
#' Additionally, a single formula may be passed instead of placing a single
#' formula in a list; e.g. `everything() ~ 1` is equivalent to
#' passing `list(everything() ~ 1)`
#'
#' @param x list of selecting formulas
#' @param type_check A predicate function that checks the elements passed on
#' the RHS of the formulas in `x=` (or the element in a named list)
#' satisfy the function.
#' @param type_check_msg When the `type_check=` fails, the string provided
#' here will be printed as the error message. When `NULL`, a generic
#' error message will be printed.
#' @param null_allowed Are `NULL` values accepted for the right hand side of
#' formulas?
#' @inheritParams .select_to_varnames
#' @keywords internal
#' @export
.formula_list_to_named_list <- function(x, data = NULL, var_info = NULL,
                                        arg_name = NULL, select_single = FALSE,
                                        type_check = NULL, type_check_msg = NULL,
                                        null_allowed = TRUE) {

  # if NULL provided, return NULL ----------------------------------------------
  if (is.null(x)) {
    return(NULL)
  }

  # converting to list if single element passed --------------------------------
  if (inherits(x, "formula")) {
    x <- list(x)
  }

  # checking the input is valid ------------------------------------------------
  .check_valid_input(x = x, arg_name = arg_name, type_check = type_check)

  # convert to a named list ----------------------------------------------------
  len_x <- length(x)
  named_list <- vector(mode = "list", length = len_x)
  for (i in seq_len(len_x)) {
    if (rlang::is_named(x[i])) {
      named_list[i] <- list(x[i])
    } else if (rlang::is_formula(x[[i]])) {
      named_list[i] <-
        .single_formula_to_list(x[[i]],
                                data = data,
                                var_info = var_info,
                                arg_name = arg_name,
                                select_single = select_single,
                                type_check = type_check,
                                type_check_msg = type_check_msg,
                                null_allowed = null_allowed
        ) |>
        list()
    } else {
      .formula_select_error(arg_name = arg_name)
    }

    .rhs_checks(
      x = named_list[i][[1]], arg_name = arg_name, type_check = type_check,
      type_check_msg = type_check_msg, null_allowed = null_allowed
    )
  }
  named_list <- purrr::flatten(named_list)

  # removing duplicates (using the last one listed if variable occurs more than once)
  rd <- function(x) {
    x <- rev(x)
    x <- !duplicated(x)
    rev(x)
  }
  tokeep <- names(named_list) |> rd()
  result <- named_list[tokeep]

  if (isTRUE(select_single) && length(result) > 1) {
    .select_single_error_msg(names(result), arg_name = arg_name)
  }
  result
}


