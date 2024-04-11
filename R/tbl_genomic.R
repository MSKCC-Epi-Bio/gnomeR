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

  # * Wide ------------

  # create data set with only alterations
  alts_only <- gene_binary %>%
    dplyr::select(-all_of(by), -sample_id)

  # create data set with all endings (including .mut if needed)
  alts_only_with_endings <- alts_only
  names(alts_only_with_endings) <- .paste_endings(names(alts_only))

  # identify types of alterations in data
  possible_alt_endings <-  c(".mut", ".Amp", ".Del", ".fus")

  which_alt_types <- purrr::map_lgl(possible_alt_endings,
                              ~any(grepl(.x, names(alts_only_with_endings))))

  names(which_alt_types) <- possible_alt_endings
  which_alt_types_no_mut <- which_alt_types[names(which_alt_types) != ".mut"]

  more_than_1_type <- (sum(which_alt_types_no_mut) > 0)

  if(wide_format == TRUE & !more_than_1_type) {

    wide_format = FALSE
    cli::cli_alert_warning("No .Amp, .Del or .fus found in your data. Ignoring `wide_format = TRUE` and returning a table of one column of summarized alterations (`wide_format = FALSE`).")

  }

# Long Table -----------------------------------------------------------


  if(wide_format == FALSE) {

    # * Order Genes for Final Table  ---------------------------------------------

    order_genes <- gene_binary %>%
      dplyr::select(-all_of(by)) %>%
      gnomeR::subset_by_frequency(t = 0) %>%
      names()

    table_data <- gene_binary %>%
      dplyr::select(all_of(by),
                    any_of(order_genes)) %>%
      dplyr::select(-"sample_id")

    # * Construct Final Table  ---------------------------------------------------

    final_table <- table_data %>%
      gtsummary::tbl_summary(by = any_of(by),...)

    if (!is.null(by)) {
      final_table <- final_table %>%
        gtsummary::add_overall()
    }

    final_table
  }

  # Wide Table -----------------------------------------------------------

  if (wide_format) {

    # * Get Gene Level Totals and Order -----------
    gene_level_binary <- summarize_by_gene(gene_binary, other_vars = any_of(by))

    gene_level_alts_only <- gene_level_binary %>%
      select(-sample_id, -any_of(by))

    order_genes <- names(gene_level_alts_only)[order(apply(gene_level_alts_only, 2, sum),
                                                   decreasing = TRUE)]

    # * Create overall n table --------------
    gene_total_table <- gene_level_alts_only %>%
      tbl_summary()


    # * Create Each Alt Type Table ------------------

    gene_binary_with_endings <- gene_binary %>%
      select(sample_id, any_of(by)) %>%
      bind_rows(alts_only_with_endings)

    # # identify types of alterations in data
    # possible_alt_endings <- c(".mut", ".Amp", ".Del", ".fus")
    #
    # # add .mut to endings to make easier
    # order_alts_explicit <- order_alts %>%
    #   .paste_endings()

    # create each table
    which_alt_types[which_alt_types == TRUE]

    tbls_list <- purrr::map(c(".mut", ".Amp", ".Del", ".fus"),
                     ~make_alt_table(alt_type = .x,
                                     data_with_endings = gene_binary_with_endings,
                                     order_genes = order_genes))

    names(tbls_list) <- c(".mut", ".Amp", ".Del", ".fus")
    tbls_list <- tbls_list %>% purrr::keep( ~ !is.null(.) )

    tbls_list$overall <- gene_total_table



  # Construct Final Table  ---------------------------------------------------
    tab_spanner_vec <- c(
      "**Overall**",
      purrr::map2(which_alt_types,
                  c("**Mutations**", "**Amplifications**", "**Deletions**", "**Fusions**"),
           ~if(.x){.y}) %>%
        unlist()
    )
        # need to figure out how to add by variables HERE
        names(data) <- c("sample_id", any_of(by),
                         .paste_endings(names(data)[2:length(names(data))]))

        data <- data %>%
          dplyr::select(any_of(by),
                        ends_with(gene))

        names(data) <- .remove_endings(names(data))

        genes_not_obs <- dplyr::setdiff(genes, names(data))[dplyr::setdiff(genes, names(data)) %in% by]

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

# Assign Class ------------------------------------------------------------


  if (wide_format) {
    class(final_table) <- c("tbl_genomic_wide", class(final_table))
  } else {
    class(final_table) <- c("tbl_genomic", class(final_table))

  }

}





# Individual Alt Type Tables ----------------------------------------------

make_alt_table <- function(alt_type = ".mut",
         data_with_endings = gene_binary_with_endings,
         order_genes) {

  any_alts <- gene_binary_with_endings  %>%
    select(ends_with(alt_type))


  if(ncol(any_alts) > 0) {
    alt_type_df <- gene_binary_with_endings  %>%
      select(sample_id, any_of(by), ends_with(alt_type))

    names(alt_type_df) <- .remove_endings(names(alt_type_df))

  # genes_not_obs <- setdiff(genes, names(data))[setdiff(genes, names(data)) %in% by]
  #
  # data <- data %>%
  # # need to fill in 0 for all unobserved alterations in hugo_symbol
  # dplyr::mutate(!!!setNames(rep(0, length(genes_not_obs)),
  # genes_not_obs))

    alt_type_df %>%
      select(any_of(order_genes)) %>%
      gtsummary::tbl_summary(by = any_of(by))
  }


}

