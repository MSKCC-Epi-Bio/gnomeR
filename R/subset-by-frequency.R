#' Subset a Binary Matrix By Alteration Frequency Threshold
#'
#' @param gene_binary A data frame with a row for each sample and column for each
#' alteration. Data frame must have a `sample_id` column and columns for each alteration
#' with values of 0, 1 or NA.
#' @param t Threshold value between 0 and 1 to subset by. Default is 10% (.1).
#' @param other_vars One or more column names (quoted or unquoted) in data to be retained
#' in resulting data frame. Default is NULL.
#' @return a data frame with a `sample_id` column and columns for
#' alterations over the given prevalence threshold of `t`.
#'
#' @export
#'
#' @examples
#' samples <- unique(gnomeR::mutations$sampleId)
#'  gene_binary <- create_gene_binary(
#'    samples = samples, mutation = mutations, cna = cna,
#'    mut_type = "somatic_only",
#'    include_silent = FALSE,
#'    specify_panel = "impact"
#'  )
#'gene_binary %>%
#'  subset_by_frequency()
#'
subset_by_frequency <- function(gene_binary, t = .1, other_vars = NULL) {


  # Checks ------------------------------------------------------------------

  # check threshold argument
  if(!(is.numeric(t) & (t >= 0 & t <= 1))) {
    cli::cli_abort("{.field t} must be a number between 0 and 1")
  }

  if (!is.data.frame(gene_binary)) {
    cli::cli_abort("{.code gene_binary} must be a data.frame")
  }

  .check_required_cols(gene_binary, "sample_id", "gene_binary")

  # Capture Other Columns to Retain -----------------------------------

  other_vars <-
    .select_to_varnames({{ other_vars }},
                        data = gene_binary,
                        arg_name = "other_vars"
    )

  # data frame of only alterations
  alt_only <- select(gene_binary, -"sample_id", -any_of(other_vars))

  # Remove all NA columns ----------------------------------------------
  all_na_alt <- apply(alt_only,  2, function(x) {
     sum(is.na(x)) == nrow(alt_only)
  })

  all_non_na_alt <- names(all_na_alt[!all_na_alt])
  alt_only <- select(alt_only, all_of(all_non_na_alt))


  # Check Numeric Class -----------------------------------------------------
  is_numeric <- apply(alt_only, 2, is.numeric)

  if(!(all(is_numeric))) {
    cli::cli_abort("All alterations in your gene binary must be numeric and only can have values of 0, 1, or NA.
                   Please coerce the following columns to numeric, or pass them to the `other_vars` argument before proceeding: {.field {names(is_numeric[!is_numeric])}}")
  }


  # Calc Frequency ----------------------------------------------------------
  counts <- apply(alt_only, 2,  function(x) {sum(x, na.rm = TRUE)})
  num_non_na <- apply(alt_only, 2, function(x) sum(!is.na(x)))

  alt_freq <- counts/num_non_na
  alts_over_thresh <- names(sort(alt_freq[alt_freq >= t], decreasing = TRUE))

  subset_binary <- select(gene_binary, "sample_id",
                          any_of(other_vars),
                          all_of(alts_over_thresh))

  return(subset_binary)

}
