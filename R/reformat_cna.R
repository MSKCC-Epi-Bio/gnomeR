
#' Reformat CNA from maf (long) version to wide version
#'
#' @param cna a cna dataframe in maf (long) format
#' @return a dataframe of reformatted CNA alteration (in wide format)
#' @export
#' @examples
#'    cna_long <- data.frame(
#'     sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
#'                  "P-0005436-T01-IM3",
#'                  "P-0001276-T01-IM3","P-0003333-T01-IM3"),
#'     Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
#'                     "HIST1H3B","KDR"),
#'     alteration = c("AMPLIFICATION","AMPLIFICATION",
#'                    "AMPLIFICATION","AMPLIFICATION","DELETION"))
#'
#' cna <- reformat_cna(cna_long)
#'
#'  cna_long <- data.frame(
#'     sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
#'                  "P-0005436-T01-IM3",
#'                  "P-0001276-T01-IM3","P-0003333-T01-IM3"),
#'     Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
#'                     "HIST1H3B","KDR"),
#'     alteration = c(2, 2, -1, 1, -2))
#'
#' cna <- reformat_cna(cna_long)
#'

reformat_cna <- function(cna) {

  cna <- cna %>%
    janitor::clean_names()


  # Check data -----------------------------------------------------------------

  missing_cols <-
    setdiff(c("sample_id", "hugo_symbol", "alteration"), names(cna))

  switch(length(missing_cols) > 0,
         cli::cli_abort("Missing columns: {.field {missing_cols}}"))

  accepted_levels <- c("NEUTRAL","LOH", "GAIN", "AMPLIFICATION", "DELETION")

  if(is.character(cna$alteration)) {
    cna <- cna %>%
      mutate(alteration = toupper(.data$alteration))

    levels_in_data <- unique(cna$alteration)

    unrecognized_coding <- setdiff(levels_in_data, accepted_levels)

    switch(length(unrecognized_coding > 0),
           cli::cli_abort("Unrecognized alteration types. Expecting {.val {accepted_levels}}"))

    # Recode Alteration Data -----------------------------------------------------
    cna <- cna %>%
      mutate(alteration =
               case_when(
                 .data$alteration == "NEUTRAL" ~ 0,
                 .data$alteration == "LOH" ~ -1,
                 .data$alteration == "GAIN" ~ 1,
                 .data$alteration == "AMPLIFICATION" ~ 2,
                 .data$alteration == "DELETION" ~ -2)
      )


  }

  hugo_symbols <- unique(cna$hugo_symbol)
  samples <- unique(cna$sample_id)



  # Create Empty Dataframe & Fill ----------------------------------------------

  # create empty data frame of correct dimensions with 0's ---
  cna_out <- as.data.frame(matrix(0L, ncol = length(samples) + 1,
                                  nrow = length(hugo_symbols)))

  colnames(cna_out) <- c("Hugo_Symbol", samples)
  cna_out[,1] <- hugo_symbols

  # fill in data frame ---
  for (i in samples) {

    alt_genes_this_sample <- cna %>%
      filter(.data$sample_id %in% i) %>%
      select("hugo_symbol") %>%
      unlist() %>%
      as.character()

    rows_to_fill <- match(alt_genes_this_sample, cna_out[, 1])
    cols_to_fill <- match(i, colnames(cna_out))

    cna_out[rows_to_fill, cols_to_fill] <-
      cna %>%
      filter(.data$sample_id %in% i) %>%
      select("alteration") %>%
      unlist() %>%
      as.numeric()
  }

  return(cna_out)

}

#' Reformat Wide CNA Data to Long
#'
#' @param wide_cna a cna dataframe in wide format (e.g. gnomeR::cna)
#' @param clean_sample_ids `TRUE` by default and function will clean
#' `sample_id` field to replace "." with "-". If `FALSE`,
#' no modification will be made to returned `sample_ids` field
#'
#' @return A long data frame of CNA events
#' @export
#'
#' @examples
#' \dontrun{
#' # NEEDS TO BE UPDATED
#' cna <- pivot_cna_longer(wide_cna = gnomeR::cna)
#' }
pivot_cna_longer <- function(wide_cna, clean_sample_ids = TRUE) {

  cna <- rename_columns(wide_cna)

  no_hugo <- select(cna, -"hugo_symbol")

  patient_sums <- apply(no_hugo, 2, sum, na.rm = TRUE)
  patient_with_sv <- patient_sums[patient_sums > 0] %>%
    names()

  cna <-  cna %>%
    select("hugo_symbol", all_of(patient_with_sv))

  cna_long <- cna %>%
    tidyr::pivot_longer(-"hugo_symbol",
                        names_to = "sample_id", values_to = "alteration")

  if(clean_sample_ids) {
    cna_long <- cna_long %>%
      mutate(sample_id = str_replace_all(.data$sample_id, fixed("."), "-"))

    cli::cli_alert_warning("Replacing all {.code .} to {.code -} in {.field sample_id} field (e.g. {.code P.0001930.T01.IM3} -> {.code P-0001930-T01-IM3}).
                   To prevent this, use argument {.code clean_sample_ids = FALSE}")
  }

  # check alteration column -----------------------------

  # Make sure hugo & alteration is character
  cna_long <- cna_long %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(str_trim(as.character(.data$alteration))))

  levels_in_data <- names(table(cna_long$alteration))

  allowed_chr_levels <- c(
    "neutral" = "0",
    "deletion" = "-2",
    "loh" = "-1.5",
    "loh" = "-1",
    "gain" = "1",
    "amplification" = "2"
  )

  all_allowed <- c(allowed_chr_levels, names(allowed_chr_levels))
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }


  suppressWarnings(
    cna_long <- cna_long %>%
      mutate(alteration = forcats::fct_recode(.data$alteration, !!!allowed_chr_levels))
  )

  cna_long


}
