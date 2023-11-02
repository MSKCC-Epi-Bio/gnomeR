#' tbl_genomic_wide
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

    # identify types of alterations in data
    any_mut <- table_data %>%
      dplyr::select(-ends_with(".Amp"),
                    -ends_with(".Del"),
                    -ends_with(".fus")) %>%
      ncol() > 0

    any_amp <- table_data %>%
      dplyr::select(ends_with(".Amp")) %>%
      ncol() > 0

    any_del <- table_data %>%
      dplyr::select(ends_with(".Del")) %>%
      ncol() > 0

    any_fus <- table_data %>%
      dplyr::select(ends_with(".fus")) %>%
      ncol() > 0

    # get names of genes
    names_genes <- gene_binary %>%
      gnomeR::summarize_by_gene() %>%
      dplyr::select(-"sample_id") %>%
      colnames() %>%
      as.vector()

    # create table of overall frequencies
    tbl1 <- gene_binary %>%
      gnomeR::summarize_by_gene() %>%
      dplyr::select(-"sample_id") %>%
      gtsummary::tbl_summary()

    # create table of mutation frequencies
    tbl_mut <- if (any_mut) {

      mut_df <- gene_binary %>%
        dplyr::select(all_of(by),
                      any_of(order_genes)) %>%
        dplyr::select(-ends_with(".Amp"),
                      -ends_with(".Del"),
                      -ends_with(".fus"))

      mut_df %>%
        dplyr::mutate(!!!setNames(rep(0, length(setdiff(names_genes, names(mut_df)))),
                                  setdiff(names_genes, names(mut_df)))) %>%
        dplyr::select(-"sample_id") %>%
        gtsummary::tbl_summary()

      # table_data %>%
      #   dplyr::select(-ends_with(".Amp"),
      #                 -ends_with(".Del"),
      #                 -ends_with(".fus")) %>%
      #   gtsummary::tbl_summary()
    }

    # create table of .Amp frequencies
    tbl_amp <- if (any_amp) {

      amp_df <- gene_binary %>%
        dplyr::select(all_of(by),
                      any_of(order_genes)) %>%
        dplyr::select(sample_id, ends_with(".Amp")) %>%
        dplyr::rename_with( ~ stringr::str_remove(., ".Amp"))

      amp_df %>%
        dplyr::mutate(!!!setNames(rep(0, length(setdiff(names_genes, names(amp_df)))),
                                  setdiff(names_genes, names(amp_df)))) %>%
        dplyr::select(-"sample_id") %>%
        gtsummary::tbl_summary()

      # table_data %>%
      #   dplyr::select(ends_with(".Amp")) %>%
      #   dplyr::rename_with( ~ stringr::str_remove(., '.Amp')) %>%
      #   gtsummary::tbl_summary()
    }

    # create table of .Del frequencies
    tbl_del <- if (any_del) {

      del_df <- gene_binary %>%
        dplyr::select(all_of(by),
                      any_of(order_genes)) %>%
        dplyr::select(sample_id, ends_with(".Del")) %>%
        dplyr::rename_with( ~ stringr::str_remove(., ".Del"))

      del_df %>%
        dplyr::mutate(!!!setNames(rep(0, length(setdiff(names_genes, names(del_df)))),
                                  setdiff(names_genes, names(del_df)))) %>%
        dplyr::select(-"sample_id") %>%
        gtsummary::tbl_summary()

      # table_data %>%
      #   dplyr::select(ends_with(".Del")) %>%
      #   dplyr::rename_with( ~ stringr::str_remove(., '.Del')) %>%
      #   gtsummary::tbl_summary()
    }

    # create table of .fus frequencies
    tbl_fus <- if (any_fus) {

      fus_df <- gene_binary %>%
        dplyr::select(all_of(by),
                      any_of(order_genes)) %>%
        dplyr::select(sample_id, ends_with(".fus")) %>%
        dplyr::rename_with( ~ stringr::str_remove(., ".fus"))

      fus_df %>%
        dplyr::mutate(!!!setNames(rep(0, length(setdiff(names_genes, names(fus_df)))),
                                  setdiff(names_genes, names(fus_df)))) %>%
        dplyr::select(-"sample_id") %>%
        gtsummary::tbl_summary()

      # table_data %>%
      #   dplyr::select(ends_with(".fus")) %>%
      #   dplyr::rename_with( ~ stringr::str_remove(., '.fus')) %>%
      #   gtsummary::tbl_summary()
    }

    # create list of tables
    tbls_list_pre <-
      list(tbl1,
           if (any_mut) {tbl_mut},
           if (any_amp) {tbl_amp},
           if (any_del) {tbl_del},
           if (any_fus) {tbl_fus})

    # drop NULL tables
    tbls_list <- tbls_list_pre %>% purrr::keep( ~ !is.null(.) )

    tab_spanner_vec <- c(
      "**Overall**",
      if (any_mut) "**Mutations**",
      if (any_amp) "**Amplifications**",
      if (any_del) "**Deletions**",
      if (any_fus) "**Fusions**"
    )

    # merge tables
    gtsummary::tbl_merge(tbls_list,
                         tab_spanner = tab_spanner_vec)

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
