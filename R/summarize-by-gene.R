#' Simplify binary matrix to one column per gene that counts any alteration type as 1
#'
#' This will reduce the number of columns in your binary matrix, and the
#' resulting data frame will have only 1 col per gene, as opposed to separate
#' columns for mutation/cna/fusion.
#'
#' @param gene_binary a 0/1 matrix of gene alterations
#'
#' @return a binary matrix with a row for each sample and one column per gene
#' @export
#'
#' @examples
#' samples <- unique(gnomeR::mutations$sampleId)[1:10]
#' gene_binary <- create_gene_binary(
#'   samples = samples, mutation = mutations, cna = cna,
#'   mut_type = "somatic_only",
#'   include_silent = FALSE,
#'   specify_panel = "IMPACT341"
#' ) %>%
#'   summarize_by_gene()
#'
summarize_by_gene <- function(gene_binary) {
  if (!is.data.frame(gene_binary)) {
    cli::cli_abort("{.code gene_binary} must be a data.frame with sample ids as {.code rownames(gene_binary)}")
  }

  if (!("sample_id" %in% names(gene_binary))) {
    gene_binary <- tibble::rownames_to_column(gene_binary, var = "sample_id")
  }


  simp_gene_binary <- gene_binary %>%
    ungroup() %>%
    tidyr::pivot_longer(-.data$sample_id) %>%
    mutate(name2 = str_remove_all(.data$name, ".Amp|.fus|.Del|.cna")) %>%
    group_by(.data$sample_id, .data$name2) %>%
    # if all NA ~ NA. If at least one non-na 1 or 0 ~ make 1 or 0
    mutate(
      num_na = sum(is.na(.data$value)),
      count = n(),
      sum = sum(.data$value, na.rm = TRUE)
    ) %>%
    mutate(simpl_val = case_when(
      count == num_na ~ NA_real_,
      sum > 1 ~ 1,
      TRUE ~ .data$sum
    )) %>%
    select("sample_id", "name2", "simpl_val") %>%
    distinct() %>%
    ungroup() %>%
    tidyr::pivot_wider(
      id_cols = .data$sample_id, names_from = .data$name2,
      values_from = .data$simpl_val
    )

  simp_gene_binary
}
