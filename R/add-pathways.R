#' Pathway Alterations
#'
#' Checks if certain oncogenic signaling pathways are altered.  Pathways were curated from [this paper](https://pubmed.ncbi.nlm.nih.gov/29625050/).  Please check for gene aliases before using.
#'
#' Sanchez-Vega, F., Mina, M., Armenia, J., Chatila, W. K., Luna, A., La, K. C., Dimitriadoy, S., Liu, D. L., Kantheti, H. S., Saghafinia, S., Chakravarty, D., Daian, F., Gao, Q., Bailey, M. H., Liang, W. W., Foltz, S. M., Shmulevich, I., Ding, L., Heins, Z., Ochoa, A., … Schultz, N. (2018). Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell, 173(2), 321–337.e10. <https://doi.org/10.1016/j.cell.2018.03.035>
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
#' @param bind_pathways a logical indicating whether pathway columns should be joined to main gene_binary or returned separately as a list item.
#' Default is TRUE and function will `bind_cols()` and return a data.frame. If FALSE a list will be returned.
#' @return a data frame: each sample is a row, columns are pathways, with values of 0/1 depending on pathway alteration status.
#' @export
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
                         bind_pathways = TRUE,
                         count_pathways_by = c("alteration", "gene")) {

  all_path <- gnomeR::pathways
  all_path_names <- names(all_path)

  # check arguments -----------------------------------------------------------

  switch(!(class(custom_pathways) %in% c("NULL", "character", "list")),
         cli::cli_abort("{.code custom_pathways} must be character vector, or list"))

  # if(!("sample_id" %in% names(gene_binary))) {
  #   gene_binary <- rownames_to_column(gene_binary, var = "sample_id")
  # }
  pathways_input <- pathways

  pathways <- pathways %>%
    purrr::when(is.null(.) ~ NULL,
                TRUE ~ match.arg(., all_path_names, several.ok = TRUE))


  not_valid <- pathways_input[!(pathways_input %in% all_path_names)]

  switch(length(not_valid) > 0,
         cli::cli_warn("Ignoring {.code {not_valid}}: not a known pathway. See {.code gnomeR::pathways}"))

  count_pathways_by <- match.arg(count_pathways_by, c("alteration", "gene"))

  all_cols <- colnames(gene_binary)
  mut_cols <- !(str_detect(all_cols, ".Amp|.Del|.fus|.cna"))


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

      # if counting by alteration but no evidence of fus/cna but fus/cna in data, warn
      # if(
      #   (any(purrr::map_lgl(custom_pathways, ~any(str_detect(.x, ".Amp|.Del|.fus|.cna")))) == FALSE) &
      #     (any(mut_cols == FALSE))
      #   ) {
      #
      #   cli::cli_warn("None of your {.code custom_pathways} have CNA (.Del/.Amp) or Fusions (.fus) but you have some in your data.
      #                 Assuming all pathway genes are mutations and counting only mutations in your data towards that pathway.
      #                 See {.code count_pathways_by} argument to change this`")
      # }

      # add mut on custom pathways when count_pathways_by == "alteration"
      if(any(purrr::map_lgl(custom_pathways, ~any(!str_detect(.x, ".Amp|.Del|.fus|.cna|.mut"))))) {
        cli::cli_inform("Assuming any gene in {.code custom_pathway} without
        suffix {.code Amp|.Del|.fus|.cna} is a mutation in that pathway. To control this behvaiour, see argument {.code count_pathways_by}")
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

  path_out <- path_out %>%
    purrr::when(
      bind_pathways ~ bind_cols(gene_binary, .),
      TRUE ~ list("gene_binary" = gene_binary,
                  "pathways" = path_out))


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
  path_alt <- gene_binary %>%
    purrr::when(
      count_pathways_by == "alteration" ~
        select(., any_of(unlist(pathway_list_item, use.names=FALSE))),
      count_pathways_by == "gene" ~
        select(., contains(unlist(pathway_list_item, use.names=FALSE)))) %>%
    mutate(sum = rowSums(., na.rm = TRUE)) %>%
    transmute('pathway_{pathway_name}' := if_else(sum >= 1, 1, 0))

  return(path_alt)
}

