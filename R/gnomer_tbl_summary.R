
#' gnomeR_tbl_summary
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
#' tb1 <- gnome_tbl_summary(data = bin.mut,cutoff = .05)
#' tb2 <- gnome_tbl_summary(data = bin.mut,gene_subset = c("KRAS", "TERT"))

gnomeR_tbl_summary <- function(data, cutoff = 0, gene_subset = NULL, by = NULL, ...){


  if(!is.data.frame(data)){
    stop("Please supply a data.frame to the data parameter.")
  }

  if(cutoff < 0 || cutoff > 1){
    stop("Please select a cutoff value between 0 and 1.")
  }

  if(!is.character(gene_subset)){
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
  }else{
    genes <- df  %>%
      select(-all_of(by)) %>%
      ungroup() %>%
      pivot_longer(-sample_id) %>%
      distinct() %>%
      group_by(name) %>%
      summarise(
        sum = sum(value, na.rm = TRUE),
        count = nrow(df) - sum(is.na(value)),
        num_na = sum(is.na(value))
      ) %>%
      mutate(perc = sum / count) %>%
      filter(perc >= cutoff) %>%
      pull(name)
  }

  # filter only those > cutoff %
  df <- df %>%
    select(sample_id, all_of(by),
           one_of(genes))


  if(is.null(by)){
    df %>%
      select(-sample_id) %>%
      tbl_summary() %>%
      bold_labels()
  }else{
    df %>%
      select(-sample_id) %>%
      tbl_summary(by = by) %>%
      add_p() %>%
      add_overall() %>%
      bold_labels() %>%
      sort_p()
  }

}
