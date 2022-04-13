#' genomic_tbl_summary
#'
#' This function will select genes based on user inputs or on frequency counts and then
#' will pass the data.frame to `gtsummary::tbl_summary()`. You can specify a `by` variable and other
#' parameters that are accepted by `gtsummary::tbl_summary()`. Note the `by` variable must be merged on to
#' onto the data before using the `by` parameter in the function.
#'
#' @param binary_matrix data.frame of genetic samples
#' @param freq_cutoff A number 0 to 1 representing the minimum percent frequency you are using to select gene's for analysis.
#' Frequencies can be calculated at gene level, or alteration level (see `freq_cutoff_by_gene`).
#' @param freq_cutoff_by_gene Logical indicating whether gene selection based on frequency % should
#' be calculated at the gene level, or the alteration level. Default is TRUE, indicating all
#' alterations (e.g. TP53/TP53.Del/TP53.Del/TP53.fus) will be aggregated at gene level to
#' determine cutoff frequencies and all alteration types for top genes will be returned.
#' @param gene_subset Specific genes you want to summarize. This takes precedent over the `freq_cutoff` parameter and will
#' return all alterations associated with that gene (including mutations, CNA, fusions).
#' @param by A variable to be passed to `gtsummary::tbl_summary()`'s by parameter
#' @param ... Additional parameters that can be passed to `gtsummary::tbl_summary()`
#'
#' @return A `tbl_summary()` object
#' @export
#'
#' @examples
#' library(gnomeR)
#' library(gtsummary)
#' library(dplyr)
#' samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' binary_matrix <- binary_matrix(samples = samples, mutation = mut, cna = cna,
#'                         mut_type = "somatic_only", snp_only = FALSE,
#'                         include_silent = FALSE,
#'                         cna_relax = TRUE, specify_panel = "no", rm_empty = FALSE)

#' tb1 <- genomic_tbl_summary(binary_matrix = binary_matrix, freq_cutoff = .05)
#' tb2 <- genomic_tbl_summary(binary_matrix = binary_matrix, gene_subset = c("KRAS", "TERT"))

genomic_tbl_summary <- function(binary_matrix,
                                freq_cutoff = 0,
                                freq_cutoff_by_gene = TRUE,
                                gene_subset = NULL,
                                by = NULL, ...){

  # check arguments & prep data ------------------------------------------------

  if(!is.data.frame(binary_matrix)){
    cli::cli_abort("{.code binary_matrix} must be a data.frame")
  }

  if(!("sample_id" %in% names(binary_matrix))) {
    binary_matrix <- rownames_to_column(binary_matrix, var = "sample_id")
  }

  if(freq_cutoff < 0 || freq_cutoff > 1){
    cli::cli_abort("Please select a {.code freq_cutoff} value between {.code 0} and {.code 1}")
  }

  # can by be string or bare variable?
  if(length(by) >1){
    stop("Please only supply 1 {.code by} variable.")
  }

  # check & assign gene subset -------------------------------------------------

  # if user passes gene_subset, we will add sufix
  gene_subset <- gene_subset %>%
    purrr::when(
      is.null(.) ~ .,

      !is.character(gene_subset) ~
        cli::cli_abort("Please supply a character vector for {.code gene_subset}"),

      length(.[(. %in% colnames(binary_matrix))]) == 0 ~
        cli::cli_abort("No genes specified in {.code gene_subset} are in your binary_matrix"),

      str_detect(., ".Amp|.Del|.fus|.cna") ~
        cli::cli_abort(
          "Detected one of the following in {.code gene_subset}: {.code '.Amp|.Del|.fus|.cna'} You may
          only pass gene names (eg. 'TP53'). To only include specific alterations, consider {.code dplyr::select(df, <alterations>)}
          before passing to {.code genomic_tbl_summary()}"),


      freq_cutoff > 0 ~
        {cli::cli_inform("You've supplied both {.code gene_subset} and {.code freq_cutoff}.
                       {.code freq_cutoff} parameter will be ignored")
          return(.)},

      # return only genes found in your data
      length(.[!(. %in% colnames(binary_matrix))]) > 0 ~
        {cli::cli_warn("The following of {.code gene_subset} are not in your data: {.code {.[!(. %in% colnames(data))]}}")
        return(.[(. %in% colnames(binary_matrix))])},
      TRUE ~ .
    )

  gene_subset <- switch(!is.null(gene_subset),
         c(gene_subset,
           paste0(gene_subset, ".Amp"),
           paste0(gene_subset, ".Del"),
           paste0(gene_subset, ".fus"),
           paste0(gene_subset, ".cna")))

  # Calc Gene Frequencies (if gene_subset is NULL) --------------------------
  gene_subset <- gene_subset %||%
    {
      binary_matrix  %>%
        select(-all_of(by)) %>%

        # if freq should be calc at gene level- simplofy matrix first
        purrr::when(
          freq_cutoff_by_gene ~ summarize_by_gene(.),
                   TRUE ~ .) %>%

        ungroup() %>%
        tidyr::pivot_longer(-.data$sample_id) %>%
        distinct() %>%
        group_by(.data$name) %>%
        summarise(
          sum = sum(.data$value, na.rm = TRUE),
          count = nrow(binary_matrix) - sum(is.na(.data$value)),
          num_na = sum(is.na(.data$value))
        ) %>%
        mutate(perc = .data$sum / .data$count) %>%
        filter(.data$perc >= freq_cutoff) %>%
        arrange(desc(.data$perc)) %>%
        pull(.data$name)
    }

  if(length(gene_subset) < 1) {
    cli::cli_abort("No genes in data set match your filter criteria (see {.code freq_cutoff})")
  }

  # if freq cutoff by gene
  gene_subset <- gene_subset %>%
    purrr::when(
    freq_cutoff_by_gene ~ c(.,
                            paste0(., ".Amp"),
                            paste0(., ".Del"),
                            paste0(., ".fus"),
                            paste0(., ".cna")),
    TRUE ~ .)


  # Select Genes and Make Table-----------------------------------------------

  table_data <-  binary_matrix %>%
    select(all_of(by),  any_of(c(gene_subset)))



  df_tbl <- table_data %>%
    gtsummary::tbl_summary(by = by) %>%
    gtsummary::bold_labels() %>%

    purrr::when(
      !is.null(by) ~
        {
          . %>%
            gtsummary::add_p() %>%
            gtsummary::add_overall()

        },
      TRUE ~ .
    )

  # should we split results by MUT/CNA/Fusion? if not at least have to fix arranging of  these
  # also its confusing when freq_by_gene is TRUE and you see low freq alts in table
  df_tbl

}
