#' Filter to final sample list
#'
#' @param df a mutation, cna, or fusion data frame to filter
#' @param samples_final a vector of sample IDs to filter to
#' @return a filtered data frame only including select samples
#' @keywords internal
#'
.filter_to_sample_list <- function(df, samples_final = NULL) {

  if (!is.null(samples_final)){
    df <- df %>%
      filter(.data$sample_id %in% samples_final)
  }
  return(df)

}

#' @param mutation Raw maf dataframe containing alteration data
#' @param include_silent Silent mutations will be removed if FALSE (default). Variant classification column is needed.
#' @return a corrected maf file or an error if problems with maf
#' @keywords internal
.check_for_silent <- function(mutation, include_silent) {

  # if include_silent FALSE, check for variant classification column -----
  if(!include_silent & !("variant_classification" %in% names(mutation))) {
    cli::cli_abort("No {.var variant_classification} column found therefore
                   silent mutations can't be removed. Please set {.code include_silent = TRUE}
                   or add a {.var variant_classification} column.")
  } else {
    return(mutation)
  }
}


#' Check for fusions in maf file
#'
#' @param mutation
#' @return a data frame if no fusions found
#' @keywords internal
.check_for_fus_in_mut <- function(mutation) {

  # Check for Fusions-  Old API used to return fusions ------
  if("variant_classification" %in% names(mutation)) {

    fusions_in_maf <- mutation %>%
      filter(.data$variant_classification %in% c("Fusion", "fusion"))

    if(nrow(fusions_in_maf) > 0) {
      cli::cli_abort("It looks like you have fusions in your mutation data frame. These need to be passed to the `fusions` argument. ")
    }
  }

  return(mutation)

}

#' Infer mutation status and assume somatic if none
#'
#' @param mutation a mutation data frame
#' @return a mutations data frame with a mutation status column
#' @keywords internal
.infer_mutation_status <- function(mutation) {

  if(!("mutation_status" %in% names(mutation))) {
    cli::cli_warn("A {.field mutation_status} column was not found. It will be assumed that
            all variants are {.val SOMATIC}, or check your data follows naming guidelines in {.code gnomer::names_df}")

    mutation <- mutation %>%
      mutate(mutation_status = "SOMATIC")
  }
  return(mutation)
}


#' Infer variant type if not present in data
#'
#' Infers variant_type from reference_allele or tumor_seq_allele data
#'
#' @param mutation data frame
#' @return a mutation data frame with a variant type column
#' @keywords internal
.infer_variant_type <- function(mutation, names_mut_dict = names_mut_dict) {

  column_names <- colnames(mutation)

  # Variant_Type ---
  if (!("variant_type" %in% column_names)) {
    if (("reference_allele" %in% column_names) & ("tumor_seq_allele2" %in% column_names)) {
      mutation %>%
        mutate(
          reference_allele = as.character(.data$reference_allele),
          tumor_seq_allele2 = as.character(.data$tumor_seq_allele2),
          variant_type = case_when(
            .data$reference_allele %in% c("A", "T", "C", "G") &
              .data$tumor_seq_allele2 %in% c("A", "T", "C", "G") ~ "SNP",
            nchar(.data$tumor_seq_allele2) < nchar(.data$reference_allele) |
              .data$tumor_seq_allele2 == "-" ~ "DEL",
            .data$reference_allele == "-" |
              nchar(.data$tumor_seq_allele2) > nchar(.data$reference_allele) ~ "INS",
            nchar(.data$reference_allele) == 2 & nchar(.data$tumor_seq_allele2) == 2 ~ "DNP",
            nchar(.data$reference_allele) == 3 & nchar(.data$tumor_seq_allele2) == 3 ~ "TNP",
            nchar(.data$reference_allele) > 3 & nchar(.data$tumor_seq_allele2) == nchar(.data$reference_allele) ~ "ONP",
            TRUE ~ "Undefined"
          )
        )


      cli::cli_warn(c("Column {.field variant_type} is missing from your data. We inferred variant types using ",
                      "{.field {dplyr::first(c(names_mut_dict['reference_allele'], 'reference_allele'), na_rm = TRUE)}} and {.field {dplyr::first(c(names_mut_dict['tumor_seq_allele_2'], 'tumor_seq_allele_2'), na_rm = TRUE)}} columns"))
    } else {
      cli::cli_abort("Column {.field variant_type} is missing from your data and {.field reference_allele} and {.field tumor_seq_allele_2}
                              columns were not available from which to infer variant type.
                              To proceed, add a column specifying {.field variant_type} (e.g. {.code mutate(<your-mutation-df>, variant_type = 'SNP')}")
    }
  }

  return(mutation)
}


#' Checks genomic input file columns to ensure column names are correct
#'
#' @param df_to_check Raw maf dataframe containing alteration data
#' @param required_cols A character specifying names of columns to check
#' @return a corrected maf file or an error if problems with maf
#' @keywords internal
#' @examples
#' gnomeR:::.clean_and_check_cols(df_to_check = gnomeR::mutations)
#'
.clean_and_check_cols <- function(df_to_check,
                                  required_cols = c("sample_id", "hugo_symbol"))  {

  df_to_check <- rename_columns(df_to_check)
  column_names <- colnames(df_to_check)

  # Check required columns & data types ------------------------------------------
  .check_required_cols(df_to_check,
                       required_cols = required_cols)

  # Make sure sample ID and hugo are character
  df_to_check <- df_to_check %>%
    mutate(across(all_of(required_cols), ~as.character(.x)))

  return(df_to_check)

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

