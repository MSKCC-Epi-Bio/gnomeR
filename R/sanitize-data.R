#' Checks MAF input columns to ensure column names are correct
#'
#' @param mutation Raw maf dataframe containing alteration data
#' @param include_silent Silent mutations will be removed if FALSE (default). Variant classification column is needed.
#' @param ... other arguments passed from create_gene_binary() (recode.aliases).
#' @return a corrected maf file or an error if problems with maf
#' @keywords internal
#' @export
#'
#' @examples
#' .clean_mutation_cols(mutation = gnomeR::mutations, include_silent = FALSE)
#'
.clean_and_check_cols <- function(df_to_check,
                                 required_cols = c("sample_id", "hugo_symbol"),
                                 data_name = NULL)  {

  mutation <- rename_columns(df_to_check)
  column_names <- colnames(df_to_check)

  # Check required columns & data types ------------------------------------------
  .check_required_cols(df_to_check, required_cols, data_name)

  # Make sure sample ID and hugo are character
  df_to_check <- df_to_check %>%
    mutate(sample_id = as.character(.data$sample_id),
           hugo_symbol = as.character(.data$hugo_symbol))

  df_to_check <- df_to_check %>%
    mutate(across(all_of(required_cols), ~as.character(.x)))

  return(df_to_check)

}

#' Checks MAF input to ensure column names are correct and renamed genes are corrected
#'
#' @param mutation Raw maf dataframe containing alteration data
#' @param include_silent Silent mutations will be removed if FALSE (default). Variant classification column is needed.
#' @param ... other arguments passed from create_gene_binary() (recode.aliases).
#' @return a corrected maf file or an error if problems with maf
#' @keywords internal
#' @export
#'
#' @examples
#' sanitize_mutation_input(mutation = gnomeR::mutations, include_silent = FALSE)
#'
sanitize_mutation_input <- function(mutation, include_silent, samples_final, ...) {

  # adding this again so this function can still be used on it's own
  mutation = clean_and_check_cols(
    df_to_check = mutation,
    required_cols = c("sample_id", "hugo_symbol"),
    data_name = "mutation"
  )

  column_names <- colnames(mutation)

  # Filter to final sample list ---------
  # * I don't think this can be NULL so maybe can remove the `if` check for NULL.
  if (!is.null(samples_final)){
    mutation <- mutation %>%
      filter(sample_id %in% samples_final)
  }

  # if include_silent FALSE, check for variant classification column -----
  if(!include_silent & !("variant_classification" %in% names(mutation))) {
    cli::cli_abort("No {.var variant_classification} column found therefore
                   silent mutations can't be removed. Please set {.code include_silent = TRUE}
                   or add a {.var variant_classification} column.")
  }



  # Check for Fusions-  Old API used to return fusions ------
  if("variant_classification" %in% column_names) {

    fusions_in_maf <- mutation %>%
      filter(.data$variant_classification %in% c("Fusion", "fusion"))

    if(nrow(fusions_in_maf) > 0) {
      cli::cli_abort("It looks like you have fusions in your mutation data frame. These need to be passed to the `fusions` argument. ")
    }
  }


  # * Check suggested columns --------

  # Mutation_Status ---
  if(!("mutation_status" %in% column_names)) {
    cli::cli_warn("A {.field mutation_status} column was not found. It will be assumed that
            all variants are {.val SOMATIC}, or check your data follows naming guidelines in {.code gnomer::names_df}")

    mutation <- mutation %>%
      mutate(mutation_status = "SOMATIC")
  }

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

      cli::cli_warn("Column {.field variant_type} is missing from your data. We inferred variant types using {.field reference_allele} and {.field tumor_seq_allele2} columns")
    } else {
      cli::cli_abort("Column {.field variant_type} is missing from your data and {.field reference_allele} and {.field tumor_seq_allele2}
                              columns were not available from which to infer variant type.
                              To proceed, add a column specifying {.field variant_type} (e.g. {.code mutate(<your-mutation-df>, variant_type = 'SNP')}")
    }
  }
  return(mutation)
}


#' Check fusion data frame to ensure columns are correct
#'
#' @param fusion a fusion data frame
#' @param ... other arguments passed from create_gene_binary()
#'
#' @return a checked data frame
#' @keywords internal
#' @export
#' @examples
#' fus <- sanitize_fusion_input(fusion = gnomeR::sv)
#'
sanitize_fusion_input <- function(fusion, samples_final)  {

  # Check required columns & data types ------------------------------------------
  # adding this again so this function can still be used on it's own
  fusion = clean_and_check_cols(
    df_to_check = fusion,
    required_cols = c("sample_id", "site_1_hugo_symbol", "site_2_hugo_symbol"),
    data_name = "fusion"
  )

  # Filter to final sample list ---------
  # * I don't think this can be NULL so maybe can remove the `if` check for NULL.
  if (!is.null(samples_final)){
    fusion <- fusion %>%
      filter(sample_id %in% samples_final)
  }

  return(fusion)
}



#' Check CNA data frame to ensure columns are correct
#'
#' @param cna a cna data frame
#' @param ... other arguments passed from create_gene_binary()
#'
#' @return a checked data frame
#' @keywords internal
#' @export
#' @examples
#'
#' cna <- sanitize_cna_input(cna = cna)
#'
sanitize_cna_input <- function(cna, samples_final, ...)  {

  # Check required columns & data types ------------------------------------------
  # adding this again so this function can still be used on it's own
  cna = clean_and_check_cols(
    df_to_check = cna,
    required_cols = c("hugo_symbol", "sample_id", "alteration"),
    data_name = "cna"
  )

  # Filter to final sample list ---------
  # * I don't think this can be NULL so maybe can remove the `if` check for NULL.
  if (!is.null(samples_final)){
    cna <- cna %>%
      filter(sample_id %in% samples_final)
  }

  # Make sure hugo & alteration is character and recode
  cna <- cna %>%
    mutate(alteration = tolower(str_trim(as.character(.data$alteration)))) %>%
    mutate(alteration = recode_cna(.data$alteration))

  return(cna)
}
