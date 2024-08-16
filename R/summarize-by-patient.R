#' Simplify binary matrix to one column per patient that counts any alteration
#' type across all samples as 1
#'
#' This will reduce the number of columns in your binary matrix, and the
#' resulting data frame will have only 1 col per gene, as opposed to separate
#' columns for mutation/cna/fusion.
#'
#' Note that if samples to the same patient were sequenced on different panels,
#' any indication of an alteration is counted as an alteration, but the absence
#' of an alteration is only defined when all sequencing panels included the gene
#' and indicated that it was not altered.
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
#'   specify_panel = "IMPACT341")
#'
#' gene_binary$patient_id = extract_patient_id(gene_binary$sample_id)
#'
#' summarize_by_patient(gene_binary)
#'
summarize_by_patient <- function(gene_binary, other_vars = NULL) {


  # Checks ------------------------------------------------------------------

  if (!is.data.frame(gene_binary)) {
    cli::cli_abort("{.code gene_binary} must be a data.frame with sample ids")
  }

  # !!! I think we should allow sample ID as input but not require it
  # .check_required_cols(
  #   gene_binary,
  #   c("sample_id"))

  .check_required_cols(
    gene_binary,
    c("patient_id"),
    add_to_message = c(i = "To extract patient IDs from IMPACT sample IDs (e.g. `P-XXXXXX-TXX-IMX`), use {.code gnomeR::extract_patient_id(data$sample_id)}"))

  # Other Vars - Capture Other Columns to Retain -----------------------------------

  other_vars <-
    .select_to_varnames({{ other_vars }},
                        data = gene_binary,
                        arg_name = "other_vars"
    )


  # Create Sample Index -----------------------------------------------------


  sample_index <- gene_binary %>%
    select("patient_id") %>%
    mutate(sample_index = paste0("samp", 1:nrow(gene_binary)))

  # data frame of only alterations

  alt_only <- as.data.frame(select(gene_binary, -"patient_id", -any_of("sample_id"), -any_of(other_vars)))

  row.names(alt_only) <- sample_index$sample_index

  # check numeric class ---------
  .abort_if_not_numeric(alt_only)

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
    filter(.data$gene %in% genes_multiple)

  if(length(genes_multiple) > 0) {
    all_bin_more <- all_bin_more %>%
      group_by(.data$gene) %>%
      summarize(across(everything(), max))
  }

  # bind together and transpose
  all_bin <- rbind(all_bin_once, all_bin_more, make.row.names = FALSE) %>%
    tibble::column_to_rownames("gene")

  all_bin <- as.data.frame(t(all_bin)) %>%
    tibble::rownames_to_column("sample_index")

  # join back to sample ID and other vars
  simp_gene_binary <- all_bin %>%
    left_join(sample_index, ., by = "sample_index") %>%
    select(-c("sample_index")) %>%
    # identify patients
    # determine number of samples per patient
    group_by(.data$patient_id) %>%
    mutate(n_samples = n()) %>%
    ungroup()

  # summarize genomic information across patients
  # separate patients w/ only 1 sample vs multiple samples to improve run time
  simp_gene_binary_pt_single <- simp_gene_binary %>%
    filter(.data$n_samples == 1)

  if (nrow(simp_gene_binary %>%
           filter(.data$n_samples > 1)) >0){

    simp_gene_binary_pt_multiple <- simp_gene_binary %>%
      filter(.data$n_samples > 1) %>%
      group_by(.data$patient_id) %>%
      summarize(across(.cols = c(everything()),
                       .fns = ~case_when(
                         # if any alteration, indicate altered
                         max(c(.x, 0), na.rm = TRUE) == 1 ~ 1,
                         # no alteration only if no NAs (no na.rm)
                         max(.x) == 0 ~ 0
                       )),
                .groups = "drop")

    simp_gene_binary_pt <- bind_rows(simp_gene_binary_pt_single,
                                     simp_gene_binary_pt_multiple) %>%
      select(-"n_samples")
    } else {
    simp_gene_binary_pt <- simp_gene_binary_pt_single %>%
      select(-"n_samples")
    }

  # !!!! Discuss this
  simp_gene_binary <- left_join(simp_gene_binary_pt,
                                gene_binary %>%
                                  select(any_of(c("patient_id", other_vars))) %>%
                                  distinct(),
                                by = "patient_id") %>%
    select("patient_id", everything())

  return(simp_gene_binary)

}
