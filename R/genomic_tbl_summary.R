
#' genomic_tbl_summary
#'
#' This function will select genes based on user inputs or on frequency counts and then
#' will pass the data.frame to `gtsummary::tbl_summary()`. You can specify a `by` variable and other
#' parameters that are accepted by `gtsummary::tbl_summary()`. Note the `by` variable must be merged on to
#' onto the data before using the `by` parameter in the function.
#'
#' @param data data.frame of genetic samples
#' @param cutoff A number 0 to 1 representing the minimum percent frequency you are using to select gene's for analysis
#' @param gene_subset Specific genes you want to summarize, take precedent over the `cutoff` parameter
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
#' bin.mut <- binary_matrix(samples = samples, mutation = mut,
#'                         mut_type = "SOMATIC", snp_only = FALSE,
#'                         include_silent = FALSE,
#'                         cna_relax = TRUE, specify_panel = "no", rm_empty = FALSE)

#' tb1 <- genomic_tbl_summary(data = bin.mut,cutoff = .05)
#' tb2 <- genomic_tbl_summary(data = bin.mut,gene_subset = c("KRAS", "TERT"))

genomic_tbl_summary <- function(data, cutoff = 0, gene_subset = NULL, by = NULL, ...){


  if(!is.data.frame(data)){
    stop("Please supply a data.frame to the data parameter.")
  }

  if(cutoff < 0 || cutoff > 1){
    stop("Please select a cutoff value between 0 and 1.")
  }

  if(!is.character(gene_subset) && !is.null(gene_subset)){
    stop("Please supply a character vector.")
  }

  if(length(by) >1){
    stop("Please only supply 1 by variable.")
  }

  df <- data %>%
    rownames_to_column(var = "sample_id")

  if(!is.null(gene_subset)){
    cli::cli_alert("Note if you specify `gene_subset` the `cutoff` parameter will be ignored.")
    genes <- gene_subset

    genes_df <- df %>%
             select(.data$sample_id, contains(genes)) %>%
             ungroup() %>%
             tidyr::pivot_longer(-.data$sample_id) %>%
             mutate(name = str_remove_all(.data$name, ".Amp|.fus|.Del")) %>%
             group_by(.data$sample_id,.data$name) %>%
             summarise(
              value = max(.data$value, na.rm = TRUE)) %>%
             ungroup() %>%
              tidyr::pivot_wider(id_cols = .data$sample_id, names_from = .data$name,
                                 values_from = .data$value)
    if(is.null(by)){
      df <- genes_df %>%
            select(-.data$sample_id)
    }else{
      df <- left_join(genes_df, select(df, .data$sample_id, all_of(by) ),
                      by = "sample_id") %>%
              select(-.data$sample_id)
    }

  }else{
    genes <- df  %>%
      select(-all_of(by)) %>%
      ungroup() %>%
      tidyr::pivot_longer(-.data$sample_id) %>%
      mutate(name = str_remove_all(.data$name, ".Amp|.fus|.Del")) %>%
      distinct() %>%
      group_by(.data$name, .data$sample_id) %>%
      summarise(value = max(.data$value, na.rm = TRUE)) %>%
      ungroup() %>%
      group_by(.data$name) %>%
      summarise(
        sum = sum(.data$value, na.rm = TRUE),
        count = nrow(df) - sum(is.na(.data$value)),
        num_na = sum(is.na(.data$value))
      ) %>%
      mutate(perc = .data$sum / .data$count) %>%
      filter(.data$perc >= cutoff) %>%
      pull(.data$name)

    if(is.null(by)){
      df <- df %>%
        select(contains(genes), .data$sample_id) %>%
        tidyr::pivot_longer(-.data$sample_id) %>%
        mutate(name = str_remove_all(.data$name, ".Amp|.fus|.Del")) %>%
        distinct() %>%
        group_by(.data$name, .data$sample_id) %>%
        summarise(value = max(.data$value, na.rm = TRUE)) %>%
        ungroup() %>%
        tidyr::pivot_wider(id_cols = .data$sample_id, names_from = .data$name,
                           values_from = .data$value) %>%
        select(-.data$sample_id)

    }else{
      df <- df %>%
            select(contains(genes), all_of(by), .data$sample_id) %>%
        tidyr::pivot_longer(c(-.data$sample_id, all_of(by))) %>%
        mutate(name = str_remove_all(.data$name, ".Amp|.fus|.Del")) %>%
        distinct() %>%
        group_by(.data$name, .data$sample_id) %>%
        summarise(value = max(.data$value, na.rm = TRUE)) %>%
        ungroup() %>%
        tidyr::pivot_wider(id_cols = c(.data$sample_id, all_of(by)),
                           names_from = .data$name,
                           values_from = .data$value) %>%
        select(-.data$sample_id)
    }

  }

  # filter only those > cutoff %


  if(is.null(by)){
    df %>%
      #select(-.data$sample_id) %>%
      gtsummary::tbl_summary() %>%
      gtsummary::bold_labels()
  }else{
    df %>%
      #select(-.data$sample_id) %>%
      gtsummary::tbl_summary(by = by) %>%
      gtsummary::add_p() %>%
      gtsummary::add_overall() %>%
      gtsummary::bold_labels() %>%
      gtsummary::sort_p()
  }

}
