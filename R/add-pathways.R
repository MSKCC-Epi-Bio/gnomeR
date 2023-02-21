#' Pathway Alterations
#'
#' Input a binary matrix of patients x genes and return a dataframe with a column per pathway
#' indicating if default or custom oncogenic signaling pathways
#' are activated in each sample. Pathways were curated
#' from [Sanchez-Vega, F et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29625050/).
#'
#' Please check for gene aliases in your dataset before using.
#'
#'
#' @param gene_binary a binary matrix from `gene_binary()`
#' @param pathways a vector of pathway names to annotate. The options are `names(gnomeR::pathways)` ("RTK/RAS", "Nrf2",
#'  "PI3K", "TGFB", "p53", "Wnt", "Myc", "Cell cycle", "Hippo", "Notch"). You can pass one pathway name, multiple pathway names, or `NULL`. By default, all
#'  pathways defined in `gnomeR::pathways` will be included. Included default pathways are alteration-specific meaning a specific type of alteration (mut/cna/fusion)
#'  is required to mark a 1 for that pathway. If you want gene-level pathways, use `summarize_by_gene()` on your binary matrix, then add
#'  pathways you want as `custom_pathways`
#' @param custom_pathways a vector of alterations to annotate as a pathway, or a list of custom pathways (see `gnomeR::pathways` as example)
#' @param count_pathways_by Must be one of the following: "alteration" (default), or "gene". This determines whether
#' any type of gene alteration should be counted towards a pathway ("gene") or only specific types of alterations should be counted towards a pathway ("alteration")
#' By default, the function assumes alteration-specific pathway annotation and all default pathways are annotated this way. If a
#' custom pathway is passed with no suffix (e.g. `custom_pathway = 'TP53'`) it will assume it is a mutation.
#' @keywords internal
#' @return a data frame: each sample is a row, columns are pathways, with values of 0/1 depending on pathway alteration status.
#' @export
#' @source Sanchez-Vega, F., Mina, M., Armenia, J., Chatila, W. K., Luna, A., La, K. C., Dimitriadoy, S., Liu, D. L., Kantheti, H. S., Saghafinia, S., Chakravarty, D., Daian, F., Gao, Q., Bailey, M. H., Liang, W. W., Foltz, S. M., Shmulevich, I., Ding, L., Heins, Z., Ochoa, A., … Schultz, N. (2018). Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell, 173(2), 321–337.e10. <https://doi.org/10.1016/j.cell.2018.03.035>
#'
#' @examples
#'
#' gene_binary <- create_gene_binary(mutation = gnomeR::mutations,
#'  cna = gnomeR::cna,
#'  fusion = gnomeR::sv)
#' pathway_df <- add_pathways(gene_binary, pathways = "Notch")
#'
add_pathways <- function(gene_binary,
                         pathways = c(names(gnomeR::pathways)),
                         custom_pathways = NULL,
                         count_pathways_by = c("alteration", "gene")) {

  all_path <- gnomeR::pathways
  all_path_names <- names(all_path)

  # check arguments -----------------------------------------------------------

  # custom pathways
  switch(!(class(custom_pathways) %in% c("NULL", "character", "list")),
         cli::cli_abort("{.code custom_pathways} must be character vector, or list"))

  .check_required_cols(gene_binary, "sample_id", "gene_binary")

  # user-specified pathways
  pathways_input <- pathways

  if(!is.null(pathways)) {
    pathways <- match.arg(pathways, all_path_names, several.ok = TRUE)
  }

  not_valid <- pathways_input[!(pathways_input %in% all_path_names)]

  switch(length(not_valid) > 0,
         cli::cli_warn("Ignoring {.code {not_valid}}: not a known pathway. See {.code gnomeR::pathways}"))

  # count pathways by
  count_pathways_by <- match.arg(count_pathways_by, c("alteration", "gene"))

  all_cols <- colnames(gene_binary)
  mut_cols <- !(str_detect(all_cols, ".Amp|.Del|.fus|.cna"))

  switch(is.null(custom_pathways) & count_pathways_by == "gene",
         cli::cli_alert_warning("Annotating the default pathways by gene may be inappropriate."))

  # custom_pathways:  can be list or vector------------------------------------
  if (!is.null(custom_pathways)) {

    # if vector, turn into a list---
    if (is.atomic(custom_pathways)) {
      empt_list <- list()
      empt_list[[1]] <- custom_pathways
      names(empt_list) <- "custom"

      custom_pathways <- empt_list

      # if is list---
    } else if (is.list(custom_pathways)) {

      #if no names, assign names
      if (is.null(names(custom_pathways))) {

        names(custom_pathways) <-
          paste0("custom_", 1:length(custom_pathways))
      }
    }


    if(count_pathways_by == "alteration") {

      # add mut on custom pathways when count_pathways_by == "alteration"
      if(any(purrr::map_lgl(custom_pathways, ~any(!str_detect(.x, ".Amp|.Del|.fus|.cna|.mut"))))) {
        cli::cli_inform("Assuming any gene in {.code custom_pathway} without
        suffix {.code .Amp|.Del|.fus|.cna} is specifically a mutation in that pathway. CNA and fusions will not be counted (e.g. TP53.Del). To control this behavior, see argument {.code count_pathways_by}")
      }


      custom_pathways <- purrr::map(custom_pathways, function(x) {
        x[!str_detect(x, ".Amp|.Del|.fus|.cna")] <-
          paste0(stringr::str_trim(
            x[!str_detect(x, ".Amp|.Del|.fus|.cna")]), ".mut")

        unique(x)

        } )

    }

    if(count_pathways_by == "gene" &
       any(purrr::map_lgl(custom_pathways, ~any(str_detect(.x, ".Amp|.Del|.fus|.cna"))))) {
      cli::cli_warn("You selected {.code count_pathways_by = 'gene'}. Ignoring {.code .Amp|.Del|.fus|.cna} in {.code custom_pathway} genes passed")
    }

    }



  # get user selected pathways (if any)
  pathways <- all_path[pathways]
  final_paths <- c(custom_pathways, pathways)

  if(count_pathways_by == "gene") {
    final_paths <- purrr::map(final_paths, ~unique(str_remove_all(.x, ".mut|.Amp|.Del|.fus|.cna|.meth")))
  }

  # prep data ------------------------------------------------------------------

  # rename mut cols -assume all non CNA/Fusion are mutations
  if (any(mut_cols)) {
    all_cols[mut_cols] <- paste0(all_cols[mut_cols], ".mut")
    colnames(gene_binary) <- all_cols
  }

  # process pathways ---------------------------------------------------------------
  path_out <- purrr::imap_dfc(final_paths, ~ .sum_alts_in_pathway(gene_binary = gene_binary,
                                                                 pathway_list_item = .x,
                                                                 pathway_name = .y,
                                                                 count_pathways_by = count_pathways_by))

  # return data  ---------------------------------------------------------------

  # remove .mut that we added
  if (any(mut_cols)) {
    all_cols[mut_cols] <- str_remove(all_cols[mut_cols], ".mut")
    colnames(gene_binary) <- all_cols
  }

  path_out <- gene_binary %>%
#    select("sample_id") %>%
    bind_cols(path_out)

  return(path_out)
}


#' Sum Alterations in a Pathway
#'
#' @param gene_binary a binary matrix (see `gene_binary()`)
#' @param pathway_list_item a named list of length 1 with pathway name as name and vector of genes as first
#' and only item in list
#' @param pathway_name name of pathway
#' @param count_pathways_by passed from `add_pathways()`
#'
#' @return a dataframe of 1 column of 0/1s indicating pathway activated yes/no
#' @keywords internal
#' @export
#'
#' @examples
#' gene_binary <- create_gene_binary(mutation = gnomeR::mutations, cna = gnomeR::cna,
#' fusion = gnomeR::sv)
#' x <- .sum_alts_in_pathway(gene_binary,
#'  pathway_list_item = gnomeR::pathways[1],
#'   pathway_name = names(gnomeR::pathways[1]),
#'     count_pathways_by = "alteration")
#'
.sum_alts_in_pathway <- function(gene_binary, pathway_list_item,
                                 pathway_name,
                                 count_pathways_by) {
  path_alt <- switch(count_pathways_by,
                     "alteration" = {
                       gene_binary %>%
                         select(any_of(unlist(pathway_list_item, use.names=FALSE)))
                     },
                     "gene" = {
                       gene_binary %>%
                         select(contains(unlist(pathway_list_item, use.names=FALSE)))
                     }) %>%
    mutate(sum = rowSums(., na.rm = TRUE)) %>%
    transmute('pathway_{pathway_name}' := if_else(sum >= 1, 1, 0))

  return(path_alt)
}

