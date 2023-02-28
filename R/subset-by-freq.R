


samples <- unique(gnomeR::mutations$sampleId)
gene_binary <- create_gene_binary(
  samples = samples, mutation = mutations, cna = cna,
  mut_type = "somatic_only",
  include_silent = FALSE,
  specify_panel = "IMPACT341"
)

subset_by_frequency <- function(gene_binary, t = .1 ) {
  if (!is.data.frame(gene_binary)) {
    cli::cli_abort("{.code gene_binary} must be a data.frame with sample ids")
  }

  .check_required_cols(gene_binary, "sample_id", "gene_binary")

  alt_only <- select(gene_binary, -sample_id)

  # Remove all NA columns ----------------------------------------------
  all_na_alt <- purrr::map_lgl(alt_only, function(x) {
    sum(is.na(x)) == nrow(alt_only)
    })
  all_non_na_alt <- names(all_na_alt[!all_na_alt])
  alt_only <- select(alt_only, all_of(all_non_na_alt))


# Check Numeric Class -----------------------------------------------------
  is_numeric <- purrr::map_lgl(alt_only, ~is.numeric(.x))
  is_numeric[!is_numeric]

  if(!(all(is_numeric))) {
    cli::cli_abort("All alterations in your gene binary must be numeric and only can have values of 0, 1, or NA.
                   Please coerce the following columns to numeric before proceeding: {.field {names(is_numeric[!is_numeric])}}")
  }


# Calc Frequency ----------------------------------------------------------
  counts <- apply(alt_only, 2, sum)
  num_na <- apply(alt_only, 2, function(x) sum(!is.na(x)))


}
