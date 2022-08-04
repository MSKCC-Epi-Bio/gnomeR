
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
#'

reformat_cna <- function(cna) {

  cna <- cna %>%
    janitor::clean_names()

  hugo_symbols <- unique(cna$hugo_symbol)
  samples <- unique(cna$sample_id)

  # Check data -----------------------------------------------------------------

  missing_cols <-
    setdiff(c("sample_id", "hugo_symbol", "alteration"), names(cna))

  switch(length(missing_cols) > 1,
         cli::cli_abort("Missing columns: {.field {missing_cols}}"))

  accepted_levels <- c("NEUTRAL","LOH", "GAIN", "AMPLIFICATION", "DELETION")


  levels_in_data <- toupper(unique(cna$alteration))

  unrecognized_coding <- setdiff(levels_in_data, accepted_levels)

  switch(length(unrecognized_coding > 0),
         cli::cli_abort("Unrecognized alteration types. Expecting {.val {accepted_levels}}"))

  # Recode Alteration Data -----------------------------------------------------
  cna <- cna %>%
    mutate(alteration =
             case_when(
               alteration == "NEUTRAL" ~ 0,
               alteration == "LOH" ~ -1,
               alteration == "GAIN" ~ 1,
               alteration == "AMPLIFICATION" ~ 2,
               alteration == "DELETION" ~ -2)
    )

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
      select(.data$hugo_symbol) %>%
      unlist() %>%
      as.character()

    rows_to_fill <- match(alt_genes_this_sample, cna_out[, 1])
    cols_to_fill <- match(i, colnames(cna_out))

    cna_out[rows_to_fill, cols_to_fill] <-
      cna %>%
      filter(.data$sample_id %in% i) %>%
      select(.data$alteration) %>%
      unlist() %>%
      as.numeric()
  }

  return(cna_out)

}
