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
    cli::cli_abort("{.code gene_binary} must be a data.frame with sample ids")
  }

  defaultW <- getOption("warn")
  options(warn = -1)

  .check_required_cols(gene_binary, "sample_id", "gene_binary")

  all_bin2 <- t(as.matrix(gene_binary))

  colnames(all_bin2) <- all_bin2[row.names(all_bin2) == "sample_id",]

  all_bin2 <- all_bin2[!row.names(all_bin2) == "sample_id", ] %>%
    as.data.frame()

  # remove endings of gene names
  all_bin2 <- all_bin2 %>%
    mutate(gene = row.names(all_bin2),
           name2 = str_remove_all(gene, ".Amp|.fus|.Del|.cna"))

  test <- all_bin2 %>%
    group_by(name2)%>%
    summarise(count = n())%>%
    filter(count > 1)

  genes_multiple <- test$name2

  # genes only with one type of event
  all_bin_once <- all_bin2 %>%
    filter(!(name2 %in% genes_multiple))%>%
    select(-gene)%>%
    mutate(across(starts_with("P"), as.numeric))

  row.names(all_bin_once) <- NULL

  # genes with more than one type of event
  all_bin_more <- all_bin2 %>%
    filter(name2 %in% genes_multiple)%>%
    mutate(across(starts_with("P"), as.numeric))%>%
    select(-gene)%>%
    group_by(name2)%>%
    summarize(across(everything(), max))

  # bind together and transpose
  all_bin <- rbind(all_bin_once, all_bin_more) %>%
    as.matrix()%>%
    t()%>%
    as.data.frame()

  colnames(all_bin) <- all_bin[row.names(all_bin) == "name2",]


  # tidy up
  all_bin <- all_bin %>%
    mutate(sample_id = row.names(.))%>%
    mutate(across(!sample_id, as.numeric))%>%
    relocate(sample_id)%>%
    filter(sample_id != "name2")

  simp_gene_binary <- all_bin

  row.names(simp_gene_binary) <- NULL


  # simp_gene_binary <- gene_binary %>%
  #   ungroup() %>%
  #   tidyr::pivot_longer(-"sample_id") %>%
  #   mutate(name2 = str_remove_all(.data$name, ".Amp|.fus|.Del|.cna")) %>%
  #   group_by(.data$sample_id, .data$name2) %>%
  #   # if all NA ~ NA. If at least one non-na 1 or 0 ~ make 1 or 0
  #   mutate(
  #     num_na = sum(is.na(.data$value)),
  #     count = n(),
  #     sum = sum(.data$value, na.rm = TRUE)
  #   )%>%
  #   mutate(simpl_val = case_when(
  #     count == num_na ~ NA_real_,
  #     sum > 1 ~ 1,
  #     TRUE ~ .data$sum
  #   ))
  #  %>%
  #   select(all_of(c("sample_id", "name2", "simpl_val"))) %>%
  #   distinct() %>%
  #   ungroup() %>%
  #   tidyr::pivot_wider(
  #     id_cols = "sample_id", names_from = "name2",
  #     values_from = "simpl_val"
  #   )

  options(warn = defaultW)

  simp_gene_binary

}










