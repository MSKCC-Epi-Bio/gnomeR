

#' IMPACT Panel Annotation of NA's
#'
#' @param gene_binary a processed gene_binary
#' @keywords internal
#' @return  a data frame iwth NAs inserted for genes not tested for given panel versions
#' @export
#'
#'
specify_impact_panels <- function(gene_binary) {

  gene_panels <- gnomeR::gene_panels

  # create data frame of sample IDs
  sample_panel_pair <- rownames(gene_binary) %>%
    as.data.frame() %>%
    stats::setNames("sample_id")

  any_impact <- sum(stringr::str_detect(sample_panel_pair$sample_id,
                                        "-IM|-IH"))

  # check if there are any impact (based on sample ID)
  switch(any_impact == 0,
         cli::cli_abort("There are no IMPACT samples recognized (based on sample_id). If you wish to annotate
                   missingness based on a panel, please pass a data.frame of sample_ids and corresponding panels."))


  # get which IMPACT panel
  sample_panel_pair <- sample_panel_pair %>%
    mutate(panel_id = case_when(
      stringr::str_detect(.data$sample_id, "-IM3") ~ "IMPACT341",
      stringr::str_detect(.data$sample_id, "-IM5") ~ "IMPACT410",
      stringr::str_detect(.data$sample_id, "-IM6") ~ "IMPACT468",
      stringr::str_detect(.data$sample_id, "-IM7") ~ "IMPACT505",
      TRUE ~ "no"
    ))

  # couldn't detect panel
  unk_impact_panel <- sample_panel_pair %>%
    filter(.data$panel_id == "no")

  if(nrow(unk_impact_panel) > 0) {
    cli::cli_alert("Couldn't infer IMPACT panel version from these sample_ids, therefore no NA panel annotation will be done for these: {unk_impact_panel$sample_id}")
  }

  return(sample_panel_pair)
}

#' Annotate Missing Gene Values According to Specific Panels
#'
#' @param sample_panel_pair a data frame of `sample_id`-`panel_id` pairs specifying panels to use for annotation of each sample
#' @param gene_binary a binary matrix of 0/1 indicating alteration yes/no for each sample
#' @keywords internal
#' @return a gene_binary annotated for missingness
#' @export

annotate_any_panel <- function(sample_panel_pair, gene_binary) {

  # if all "no", leave function
  switch(all(sample_panel_pair$panel_id == "no"),
         return(gene_binary))

  sample_panel_pair_nest <- sample_panel_pair %>%
    group_by(.data$panel_id) %>%
    summarise(samples_in_panel = list(.data$sample_id))

  # pull genes for given panels
  panels_needed <- unique(sample_panel_pair_nest$panel_id)

  # has sample IDs and genes for each panel
  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    left_join(gnomeR::gene_panels, by = c("panel_id" = "gene_panel")) %>%
    select(-"entrez_ids_in_panel")

  user_data_genes <- gsub(".fus|.Del|.Amp|.cna", "", colnames(gene_binary))

  sample_panel_pair_nest <- sample_panel_pair_nest %>%
    mutate(na_genes_raw = purrr::map(.data$genes_in_panel,
                                     ~unique(setdiff(user_data_genes, .x)))) %>%
    mutate(na_genes = purrr::map(.data$na_genes_raw,
                                 ~c(
                                   .x,
                                   paste0(.x, ".fus"),
                                   paste0(.x, ".Del"),
                                   paste0(.x, ".Amp"),
                                   paste0(.x, ".cna")
                                 )))


  annotated_data <- purrr::pmap_df(sample_panel_pair_nest,
                                   annotate_specific_panel,
                                   gene_binary = gene_binary)

  return(annotated_data)
}


#' Utility function  to insert NA's According to Panel
#'
#' @param gene_binary a processed binary matrix
#' @param panel_id name of gene panel
#' @param samples_in_panel samples to be annotated for each panel
#' @param na_genes genes to make NA
#' @param ... other args passed
#' @keywords internal
#' @return an annotated data frame
#' @export
#'
#'
annotate_specific_panel <- function(gene_binary,
                                    panel_id,
                                    samples_in_panel,
                                    na_genes, ...) {

  mut_sub <- gene_binary[samples_in_panel, ]
  mut_sub[,stats::na.omit(match(na_genes, colnames(mut_sub)))] <- NA

  return(mut_sub)

}


# Check which panels gene is on----------------------

#' provide a list of impact panels a provided gene is found within
#'
#' @param hugo_symbol a vector of hugo symbols
#'
#' @return a data frame with hugo symbols and the IMPACT panels on which they are
#' included
#' @keywords internal
#' @examples
#'
#' hugos <- unique(gnomeR::mutations$hugoGeneSymbol)[1:10]
#'
#' which_impact_panel(hugos)
#'
#' @export

which_impact_panel <- function(hugo_symbol) {

  # get table of gene aliases (internal data)
  alias_table <- gnomeR::impact_alias_table %>%
    dplyr::select("hugo_symbol", "alias")

  # recode all genes to most common alias
  hugo_symbol <- purrr::map_chr(hugo_symbol,
                                ~resolve_alias(gene_to_check = .x,
                                               alias_table = alias_table))

  gene_panels <- gnomeR::gene_panels %>%
    filter(.data$gene_panel %in% c("IMPACT341", "IMPACT410", "IMPACT468", "IMPACT505")) %>%
    select("gene_panel", "genes_in_panel") %>%
    tidyr::unnest(cols = c("genes_in_panel"))


  impact_results <- gene_panels %>%
    filter(.data$genes_in_panel %in% hugo_symbol) %>%
    distinct() %>%
    mutate(fill = "yes") %>%
    tidyr::pivot_wider(
      names_from = "gene_panel",
      values_from = "fill",
      values_fill = "no")

  impact_results

}
