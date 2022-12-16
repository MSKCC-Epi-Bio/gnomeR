
#' Pivot CNA from maf (long) version to wide version
#'
#' @param cna a cna dataframe in maf (long) format
#' @return a dataframe of reformatted CNA alteration (in wide format)
#' @export
#' @examples
#' cna_long <- data.frame(
#'     sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
#'                  "P-0005436-T01-IM3",
#'                  "P-0001276-T01-IM3","P-0003333-T01-IM3"),
#'     Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
#'                     "HIST1H3B","KDR"),
#'     alteration = c("AMPLIFICATION","AMPLIFICATION",
#'                    "AMPLIFICATION","AMPLIFICATION","DELETION"))
#'
#' cna <- pivot_cna_wider(cna_long)
#'
#'  cna_long <- data.frame(
#'     sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
#'                  "P-0005436-T01-IM3",
#'                  "P-0001276-T01-IM3","P-0003333-T01-IM3"),
#'     Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
#'                     "HIST1H3B","KDR"),
#'     alteration = c(2, 2, -1, 1, -2))
#'
#' cna <- pivot_cna_wider(cna_long)
#'

pivot_cna_wider <- function(cna) {

  cna <- cna %>%
    janitor::clean_names()


  # Check data -----------------------------------------------------------------

  missing_cols <-
    setdiff(c("sample_id", "hugo_symbol", "alteration"), names(cna))

  switch(length(missing_cols) > 0,
         cli::cli_abort("Missing columns: {.field {missing_cols}}"))

  cna <- cna %>%
    mutate(alteration = as.character(.data$alteration)) %>%
    mutate(alteration = toupper(.data$alteration))

  accepted_levels <- c("NEUTRAL","LOH", "GAIN", "AMPLIFICATION", "DELETION",
                       "-2", "-1", "0", "1", "2")

  levels_in_data <- unique(cna$alteration)

  not_allowed <- stats::na.omit(setdiff(levels_in_data, accepted_levels))

  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {accepted_levels}}"))
  }

  if(any(cna$alteration %in% c("NEUTRAL","LOH", "GAIN", "AMPLIFICATION", "DELETION"))) {

    # Recode Alteration Data -----------------------------------------------------
    cna <- suppressWarnings(cna %>%
      mutate(alteration =
               case_when(
                 .data$alteration == "NEUTRAL" ~ 0,
                 .data$alteration == "LOH" ~ -1,
                 .data$alteration == "GAIN" ~ 1,
                 .data$alteration == "AMPLIFICATION" ~ 2,
                 .data$alteration == "DELETION" ~ -2,
                 TRUE ~ as.numeric(.data$alteration))))

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
#' Takes a numeric vector of CNA data in wide format where each column is a sample
#' and each row is a hugo symbol. Function will return a long format CNA data set
#' of just events (neutral/diploid instances are filtered out) and will recode events from
#' numeric to descriptive (-2/-1/-1.5 is deletion, 2/1 is amplification).
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
#' cna <- pivot_cna_longer(wide_cna = gnomeR::cna_wide)

pivot_cna_longer <- function(wide_cna, clean_sample_ids = TRUE) {

  cna <- rename_columns(wide_cna)

  # Check data -----------------------------------------------------------------

  missing_cols <-
    setdiff(c("hugo_symbol"), names(cna))

  switch(length(missing_cols) > 0,
         cli::cli_abort("Missing columns: {.field {missing_cols}}"))

  no_hugo <- select(cna, -"hugo_symbol")

  switch(any(purrr::map_lgl(no_hugo, ~!is.numeric(.x))),
         cli::cli_abort(c("All CNA columns must be numeric. Do you need to convert to numeric? Eg. ",
         "{.code mutate(data, across(.cols = c(everything(), -Hugo_Symbol), ~ as.numeric(.x)))} ?"))
  )

  # Remove Patients with no CNA ------------------------------------------------


  patient_sums <- apply(no_hugo, 2, function(x) sum(abs(x), na.rm = TRUE))

  patient_with_cna <- patient_sums[patient_sums > 0] %>%
    names()

  if(length(patient_with_cna) < 1) {
    cli::cli_abort("There are no CNA events in the data set.")
  }

  # Pivot Longer ---------------------------------------------------------------
  cna <-  cna %>%
    select("hugo_symbol", all_of(patient_with_cna))

  cna_long <- cna %>%
    tidyr::pivot_longer(-"hugo_symbol",
                        names_to = "sample_id", values_to = "alteration")

  # remove neutral events (pre-filtered above for speed but may not need this)
  cna_long <- cna_long %>%
    filter(.data$alteration != 0)

  # clean Names  ---------------------------------------------------------------
  if(clean_sample_ids) {
    cna_long <- cna_long %>%
      mutate(sample_id = str_replace_all(.data$sample_id, fixed("."), "-"))

    cli::cli_alert_warning("Replacing all {.code .} to {.code -} in {.field sample_id} field (e.g. {.code P.0001930.T01.IM3} -> {.code P-0001930-T01-IM3}).
                   To prevent this, use argument {.code clean_sample_ids = FALSE}")
  }

  # Recode Alteration ---------------------------------------------------------

  # Make sure hugo & alteration is character
  cna_long <- cna_long %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(str_trim(as.character(.data$alteration))))

  # recode alterations
  cna_long <- cna_long %>%
    mutate(alteration = recode_cna(.data$alteration))

  return(cna_long)


}
