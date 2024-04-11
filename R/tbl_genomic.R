#' tbl_genomic
#'
#' This function will select genes based on user inputs or on frequency counts
#' and then will pass the data.frame to `gtsummary::tbl_summary()`. You can
#' specify a `by` variable and other parameters that are accepted by
#' `gtsummary::tbl_summary()`. Note the `by` variable must be merged on to onto
#' the data before using the `by` parameter in the function.
#'
#' @param gene_binary data.frame of genetic samples
#' @param freq_cutoff deprecated
#' @param freq_cutoff_by_gene deprecated
#' @param gene_subset deprecated
#' @param wide_format Specifies whether to stratify tbl_genomic by alteration
#'   type such that the resulting table will include one column per alteration
#'   type and an overall summary column. Default is `wide_format = FALSE`.
#' @param by A variable to be passed to `gtsummary::tbl_summary()`'s by
#'   parameter
#' @param ... Additional parameters that can be passed to
#'   `gtsummary::tbl_summary()`. To access the additional parameters you need to
#'   load `gtsummary`.
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
                        wide_format = FALSE,
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


  # * Wide format by alteration type -----

  if (wide_format) {

    final_table <- .create_wide_format(
      gene_binary_df = gene_binary,
      order_genes_df = order_genes,
      table_data_df = table_data
    )

  }

  # Construct Final Table  ---------------------------------------------------

  else {final_table <- table_data %>%
    gtsummary::tbl_summary(by = any_of(by),...)

  if (!is.null(by)) {
    final_table <- final_table %>%
      gtsummary::add_overall()
  }

  final_table}

}

#' * Wide format by alteration type -----
#'
#' This helper function will stratify tbl_genomic by alteration type such that
#' the resulting table will include one column per alteration type and an
#' overall summary column.
#'
#' @param gene_binary_df data.frame of genetic samples
#' @param order_genes_df data.frame
#' @param table_data_df data.frame
#'
#' @return wide format 'tbl_genomic' by alteration type
#' @keywords internal
#' @export
#'
.create_wide_format <- function(gene_binary_df, order_genes_df, table_data_df){

  genes <- order_genes_df[!order_genes_df %in% by]%>%
    .remove_endings()%>%
    unique()

  # identify types of alterations in data

  gene_endings <- c(".mut", ".Amp", ".Del", ".fus")

  # add .mut to endings to make easier

  gb_alt_all <- order_genes_df %>%
    .paste_endings()

  any_alt_types <- purrr::map(gene_endings,
                              ~any(grepl(.x, gb_alt_all)))

  names(any_alt_types) <- gene_endings

  # create table of overall frequencies

  tbl1 <- gene_binary_df %>%
    tidyr::pivot_longer(!sample_id)%>%
    dplyr::mutate(name = .remove_endings(name))%>%
    dplyr::group_by(sample_id, name)%>%
    dplyr::slice(which.max(value))%>%
    dplyr::ungroup()%>%
    dplyr::select(-"sample_id")


  tbl2 <- tbl1 %>%
    split(tbl1$name)

  tbl_overall <- purrr::map(1:length(names(tbl2)), function(x){
    tbl2[[x]] %>%
      dplyr::select(-"name")%>%
      stats::setNames(names(tbl2)[[x]])
  })%>%
    do.call(cbind, .)%>%
    gtsummary::tbl_summary(by = any_of(by))

  # create table of mutation frequencies

  tbls_alt_types <- purrr::map2(
    any_alt_types,
    gene_endings,
    function(yes_exist, gene) {

      if(yes_exist[[1]]){
        data <- table_data

        # need to figure out how to add by variables HERE
        names(data) <- c("sample_id", any_of(by),
                         .paste_endings(names(data)[2:length(names(data))]))

        data <- data %>%
          dplyr::select(any_of(by),
                        ends_with(gene))

        names(data) <- .remove_endings(names(data))

        genes_not_obs <- dplyr::setdiff(genes, names(data))[dplyr::setdiff(genes, names(data)) %in% by]

        data <- data %>%
          # need to fill in 0 for all unobserved alterations in hugo_symbol
          dplyr::mutate(!!!stats::setNames(rep(0, length(genes_not_obs)),
                                           genes_not_obs))

        data %>%
          gtsummary::tbl_summary(by = any_of(by))
      } else {
        NULL
      }

    }
  )

  # create list of tables
  tbls_list_wide_pre <-
    append(list(tbl_overall),
           purrr::map2(any_alt_types,
                       tbls_alt_types, function(x, y){
                         if (x) {y}
                       }))

  # drop NULL tables
  tbls_list_wide <- tbls_list_wide_pre %>% purrr::keep( ~ !is.null(.) )


  tab_spanner_vec <- c(
    "**Overall**",
    purrr::map2(any_alt_types, c("**Mutations**", "**Amplifications**",
                                 "**Deletions**", "**Fusions**"),
                ~if(.x){.y}) %>%
      unlist()
  )

  # merge tables
  final_table_wide <- gtsummary::tbl_merge(tbls_list_wide,
                                           tab_spanner = tab_spanner_vec)

  return(final_table_wide)

}

