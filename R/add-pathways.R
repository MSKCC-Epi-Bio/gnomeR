#' Pathway Alterations
#'
#' Input a binary matrix of patients x alterations and return a dataframe with a column per pathway
#' indicating if default or custom oncogenic signaling pathways
#' are activated in each sample. Default package pathways were sourced
#' from [Sanchez-Vega, F et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29625050/).
#'
#' Please check for gene aliases in your data set before using.
#'
#'
#' @description
#' Please note that only `sample_id column`, and columns with .Amp, .Del, .fus or no suffix are accepted.
#' Any gene column with no suffix will be assumed to be a mutation.
#'
#' @param gene_binary a binary matrix from `create_gene_binary()`
#' @param pathways a vector of pre-coded pathways to annotate. The options are `names(gnomeR::pathways)` ("RTK/RAS", "Nrf2",
#'  "PI3K", "TGFB", "p53", "Wnt", "Myc", "Cell cycle", "Hippo", "Notch"). You can pass multiple pathway names, or `NULL`. By default, all
#'  pathways defined in `gnomeR::pathways` will be included. Included default pathways are alteration-specific, meaning a specific type of alteration (mut/cna/fusion)
#'  is required to mark a 1 for that pathway.
#' @param custom_pathways a vector of alterations to annotate as a single pathway, or a list of custom pathways (see `gnomeR::pathways` as example).
#' You must specify the alteration type for each gene using `.mut`, `.Amp`, `.Del` suffix, e.g. `c("TP53.mut", "CDKN2A.Amp")`. If you wish to count any type of
#' alteration on that gene towards the pathway you can use the `.any` suffix (e.g. `c("TP53.any")`).
#' @param other_vars One or more column names (quoted or unquoted) in data to be retained
#' in resulting data frame. Default is NULL.
#' @param count_pathways_by deprecated
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
                         other_vars = NULL,
                         count_pathways_by = deprecated()) {

  # Check arguments -----------------------------------------------------------

  all_path <- gnomeR::pathways
  all_path_names <- names(all_path)
  .check_required_cols(gene_binary, "sample_id")

  # * Deprecated Arguments (will remove this in the future) ----

  if (lifecycle::is_present(count_pathways_by)) {
    lifecycle::deprecate_stop(when = "1.3.0",
                              what = "add_pathways(count_pathways_by)",
                              details = c(
                                i = c("All pathways are now counted on the specific alteration level, not the gene level and all columns in your data with no suffix (e.g. `data$TP53`) are assumed to be mutations.
                                You must explicitely specify in your custom pathway what types of alterations count towards pathway e.g. `custom_pathway = c('TP53.mut', 'APC.Del').
                                      If you would like to count any type of alteration on a gene towards a pathway, use the suffix `.any`, e.g. `custom_pathway = c('TP53.any', 'APC.any')")))
  }



  # * Default pathways ----
  pathways_input <- pathways

  if(!is.null(pathways)) {
    pathways <- match.arg(pathways, all_path_names, several.ok = TRUE)
  }

  not_valid <- pathways_input[!(pathways_input %in% all_path_names)]

  switch(length(not_valid) > 0,
         cli::cli_warn("Ignoring {.code {not_valid}}: not a known pathway. See {.code gnomeR::pathways}"))


  # Process Custom Pathways ----

  switch(!(class(custom_pathways) %in% c("NULL", "character", "list")),
         cli::cli_abort("{.code custom_pathways} must be character vector, or list"))


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

    if(any(purrr::map_lgl(custom_pathways, ~any(!str_detect(.x, ".Amp|.Del|.fus|.mut|.any"))))) {
      cli::cli_abort(c("All alterations specified in {.code custom_pathway} must have one ",
                      "of the following suffixes specified: {.code .Amp}, {.code .Del}, {.code .fus}, {.code .mut}, {.code .any}"))
    }


    # process GENE.any suffix in custom pathways
    custom_pathways <- purrr::map(custom_pathways, function(x) {
      gene_all <- x[str_detect(x, ".any")]

      if(length(gene_all) > 0) {
        gene_all <-  str_remove(gene_all, ".any")

        gene_add <- purrr::map(gene_all, function(x) {
          paste0(x, c(".Amp", ".Del", ".fus", ".mut"))
        }) %>% unlist()

        x <- c(x[!str_detect(x, ".any")], gene_add)
      }

      unique(x)

      } )

    }


  # get user selected pathways (if any)
  pathways <- all_path[pathways]
  final_paths <- c(custom_pathways, pathways)

  # Prep data ------------------------------------------------------------------
  # * Other Vars - Capture Other Columns to Retain ----------------

  other_vars <-
    .select_to_varnames({{ other_vars }},
                        data = gene_binary,
                        arg_name = "other_vars"
    )

  # data frame of only alterations
  alt_only <- select(gene_binary, -"sample_id", -any_of(other_vars))


  all_cols <- colnames(alt_only)
  mut_cols <- !(str_detect(all_cols, ".Amp|.Del|.fus"))

  # * Rename .mut columns (assume all non CNA/Fusion are mutations) -----

  if (any(mut_cols)) {

    # in case any columns already had .mut suffix
    all_cols[mut_cols] <- str_remove(all_cols[mut_cols], ".mut")

    # now add .mut
    all_cols[mut_cols] <- paste0(all_cols[mut_cols], ".mut")

    colnames(alt_only) <- all_cols
  }

  # process pathways ---------------------------------------------------------------
  path_out <- purrr::imap_dfc(final_paths, ~ .sum_alts_in_pathway(gene_binary = alt_only,
                                                                 pathway_list_item = .x,
                                                                 pathway_name = .y))

  # return data  ---------------------------------------------------------------

  # # remove .mut that we added
  # if (any(mut_cols)) {
  #   all_cols[mut_cols] <- str_remove(all_cols[mut_cols], ".mut")
  #   colnames(alt_only) <- all_cols
  # }

  path_out <- gene_binary %>%
    bind_cols(path_out)

  return(path_out)
}


#' Sum Alterations in a Pathway
#'
#' @param gene_binary a binary matrix (see `gene_binary()`)
#' @param pathway_list_item a named list of length 1 with pathway name as name and vector of genes as first
#' and only item in list
#' @param pathway_name name of pathway
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
  path_alt <- gene_binary %>%
    select(any_of(unlist(pathway_list_item, use.names=FALSE))) %>%
    mutate(sum = rowSums(., na.rm = TRUE)) %>%
    transmute('pathway_{pathway_name}' := if_else(sum >= 1, 1, 0))

  return(path_alt)
}

