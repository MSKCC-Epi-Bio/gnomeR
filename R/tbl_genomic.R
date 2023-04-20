#' tbl_genomic
#'
#' This function will select genes based on user inputs or on frequency counts and then
#' will pass the data.frame to `gtsummary::tbl_summary()`. You can specify a `by` variable and other
#' parameters that are accepted by `gtsummary::tbl_summary()`. Note the `by` variable must be merged on to
#' onto the data before using the `by` parameter in the function.
#'
#' @param gene_binary data.frame of genetic samples
#' @param freq_cutoff deprecated
#' @param freq_cutoff_by_gene deprecated
#' @param gene_subset deprecated
#' @param by A variable to be passed to `gtsummary::tbl_summary()`'s by parameter
#' @param ... Additional parameters that can be passed to `gtsummary::tbl_summary()`. To access the additional parameters you need to load `gtsummary`.
#' @return A `tbl_summary()` object
#' @export
#'
#' @examples
#'
#' samples <- unique(mutations$sampleId)[1:10]
#'
#' gene_binary <- create_gene_binary(
#'   samples = samples,
#'   mutation = gnomeR::mutations,
#'   cna = gnomeR::cna,
#'   mut_type = "somatic_only", snp_only = FALSE,
#'   specify_panel = "no"
#' )
#'
#' tbl1 <- tbl_genomic(gene_binary)
#'
#' # Example wth `by` variable
#'
#' gene_binary$sex <- sample( c("M", "F"), size = nrow(gene_binary), replace = TRUE)
#'
#' tbl2 <- tbl_genomic(
#'   gene_binary = gene_binary,
#'   by = sex
#' ) %>%
#' gtsummary::add_p() %>%
#' gtsummary::add_q()
#'
tbl_genomic <- function(gene_binary,
                        by = NULL,
                        freq_cutoff = deprecated(),
                        freq_cutoff_by_gene = deprecated(),
                        gene_subset = deprecated(),
                        ...) {

  # Check arguments & prep data ------------------------------------------------

  if (!inherits(gene_binary, "data.frame")) {
    stop("`gene_binary=` argument must be a tibble or data frame.", call. = FALSE)
  }

  .check_required_cols(gene_binary, "sample_id", "gene_binary")

  if("sample_id" %in% names(gene_binary)) {
    if(any(table(gene_binary$sample_id) > 1)) {
      cli::cli_abort("Duplicate `sample_ids` found in `gene_binary`. Samples IDs should be unique.")
    }

  }
  # * Deprecated Arguments (will remove this in the future) ----

  if (lifecycle::is_present(freq_cutoff)) {
    lifecycle::deprecate_stop(when = "1.3.0",
                              what = "tbl_genomic(freq_cutoff)",
                              details = c(
                                i = c("Please pre-select gene columns to summarize before passing to `tbl_genomic()`.",
                                      "Use `gnomeR::subset_by_frequency(t)` to easily subset by a gene prevalance threshold.")))
  }

  if (lifecycle::is_present(freq_cutoff_by_gene)) {
    lifecycle::deprecate_stop(when = "1.3.0",
                              what = "tbl_genomic(freq_cutoff_by_gene)",
                              details = c(
                                i = c("Please pre-select gene columns to summarize before passing to `tbl_genomic()`.",
                                      "Use `gnomeR::subset_by_frequency(t)` to easily subset by a gene prevalance threshold.")))
  }

  if (lifecycle::is_present(gene_subset)) {
    lifecycle::deprecate_stop(when = "1.3.0",
                              what = "tbl_genomic(gene_subset)",
                              details = c(
                                i = c("Please pre-select gene columns to summarize before passing to `tbl_genomic()`.",
                                      "Use `gnomeR::subset_by_frequency(t)` to easily subset by a gene prevalance threshold.")))
  }

  # * Other Args -----

  by <-
    .select_to_varnames({{ by }},
      data = gene_binary,
      arg_name = "by", select_single = TRUE
    )

  # Order Genes for Final Table  ---------------------------------------------

  order_genes <- gene_binary %>%
    dplyr::select(-all_of(by)) %>%
    gnomeR::subset_by_frequency(t = 0) %>%
    names()

  table_data <- gene_binary %>%
    dplyr::select(all_of(by),
                  any_of(order_genes)) %>%
    dplyr::select(-"sample_id")

  # Construct Final Table  ---------------------------------------------------

  final_table <- table_data %>%
    gtsummary::tbl_summary(by = any_of(by),...)

  if (!is.null(by)) {
    final_table <- final_table %>%
      gtsummary::add_overall()
  }

  final_table

}
