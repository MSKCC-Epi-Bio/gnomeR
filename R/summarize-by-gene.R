#' Simplify binary matrix to one column per gene that counts any alteration type as 1
#'
#' This will reduce the number of columns in your binary matrix, and the
#' resulting data frame will have only 1 col per gene, as opposed to separate
#' columns for mutation/cna/fusion.
#'
#' @param gene_binary a 0/1 matrix of gene alterations
#' @param other_vars One or more column names (quoted or unquoted) in data to be retained
#' in resulting data frame. Default is NULL.
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
summarize_by_gene <- function(gene_binary, other_vars = NULL) {


  # Checks ------------------------------------------------------------------

  if (!is.data.frame(gene_binary)) {
    cli::cli_abort("{.code gene_binary} must be a data.frame with sample ids")
  }

  .check_required_cols(gene_binary, "sample_id", "gene_binary")

  # check for repeat samples
  if(any(table(gene_binary$sample_id) > 1)) {
    cli::cli_abort("Your {.field gene_binary} must have unique samples in {.code sample_id} column")
  }

  # Capture Other Columns to Retain -----------------------------------

  other_vars <-
    .select_to_varnames({{ other_vars }},
                        data = gene_binary,
                        arg_name = "other_vars"
    )


  # Create Sample Index -----------------------------------------------------

  sample_index <- gene_binary %>%
    select("sample_id") %>%
    mutate(sample_index = paste0("samp", 1:nrow(gene_binary)))

<<<<<<< Updated upstream
  alt_only <- as.matrix(select(gene_binary, -"sample_id"))
  rownames(alt_only) <- sample_index$sample_index
=======
  # data frame of only alterations
  alt_only <- select(gene_binary, -"sample_id", -any_of(other_vars)) %>%
    as.matrix()

  row.names(alt_only) <- sample_index$sample_index
>>>>>>> Stashed changes

  # check numeric class ---------
  is_numeric <- apply(alt_only, 2, is.numeric)

  if(!(all(is_numeric))) {
    cli::cli_abort("All alterations in your gene binary must be numeric and only can have values of 0, 1, or NA.
                   Please coerce the following columns to numeric before proceeding: {.field {names(is_numeric[!is_numeric])}}")
  }

  # Transpose ---------------------------------------------------------------

  transp_alt_only <- as.data.frame(t(alt_only))

  # remove endings of gene names
  transp_alt_only <- transp_alt_only %>%
    mutate(gene = str_remove_all(row.names(.),
                                 ".Amp|.fus|.Del|.cna"))

  # check for genes that have more than one alt type
  gene_tab <- table(transp_alt_only$gene)

  genes_multiple <-  names(gene_tab[which(gene_tab > 1)])
  genes_single <- names(gene_tab[which(gene_tab == 1)])

  # genes with one type of event
  all_bin_once <- transp_alt_only %>%
    filter(.data$gene %in% genes_single)

  # genes with more than one type of event
  all_bin_more <- transp_alt_only %>%
    filter(.data$gene %in% genes_multiple) %>%
    group_by(.data$gene) %>%
    summarize(across(everything(), max))

  # bind together and transpose
  all_bin <- rbind(all_bin_once, all_bin_more, make.row.names = FALSE) %>%
    tibble::column_to_rownames("gene")

  all_bin <- as.data.frame(t(all_bin)) %>%
    tibble::rownames_to_column("sample_index")

  # join back to sample ID and other vars
  simp_gene_binary <- all_bin %>%
    left_join(sample_index, ., by = "sample_index") %>%
    select(-c("sample_index")) %>%
    as.data.frame()%>%
    left_join(gene_binary %>%
                select(any_of(c("sample_id", other_vars))))%>%
    suppressMessages()

  simp_gene_binary

}










