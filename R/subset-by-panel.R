#' Subset a Binary Matrix By Genes Available on Specified Panel
#'
#' @param gene_binary A data frame with a row for each sample and column for each
#' alteration. Data frame must have a `sample_id` column and columns for each alteration
#' with values of 0, 1 or NA.
#' @param panel_id A character string or vector of the specified panel to subset the genes (see `gnomeR::gene_panels` for available panels)
#' @param other_vars One or more column names (quoted or unquoted) in data to be retained
#' in resulting data frame. Default is NULL.
#' @return a data frame with a `sample_id` column and columns for
#' alterations on genes that were sequenced on the specified panel.
#' @author Jessica Lavery
#'
#' @export
#'
#' @examples
#' samples <- unique(gnomeR::mutations$sampleId)
#' gene_binary <- create_gene_binary(
#'   samples = samples, mutation = mutations, cna = cna,
#'   mut_type = "somatic_only",
#'   include_silent = FALSE,
#'   specify_panel = "impact"
#' )
#' subset_by_panel(gene_binary = gene_binary, panel_id = "IMPACT468")
#'
# get names of genes that were on the panel but that were unaltered in this sample
#' p_genes <- tidyr::unnest(gnomeR::gene_panels, cols = c("genes_in_panel"))
#' p_genes <- p_genes[p_genes$gene_panel == 'IMPACT300', ]
#' setdiff(p_genes$genes_in_panel, names(gene_binary))

subset_by_panel <- function(gene_binary, panel_id = NULL, other_vars = NULL){

  gene_panels <- gnomeR::gene_panels

  # Checks ------------------------------------------------------------------
  # check inputs
  if (!is.data.frame(gene_binary)) {
    cli::cli_abort("{.code gene_binary} must be a data.frame")
  }

  if(is.null(panel_id)) {
    cli::cli_abort("{.code panel_id} must not be NULL")
  }

  .check_required_cols(gene_binary, "sample_id")

  if (!(panel_id %in% c(gene_panels$gene_panel))){
    cli::cli_abort("The panel {panel_id} is not an available panel. See `gnomeR::gene_panels()` for the names of available panels.")
  }

  # Other Vars - Capture Other Columns to Retain -----------------------------------
  other_vars <-
    .select_to_varnames({{ other_vars }},
                        data = gene_binary,
                        arg_name = "other_vars"
    )

  # data frame of only alterations
  alt_only <- select(gene_binary, -"sample_id", -any_of(other_vars))

  # Remove all NA columns ----------------------------------------------
  all_na_alt <- apply(alt_only,  2, function(x) {
    sum(is.na(x)) == nrow(alt_only)
  })

  all_non_na_alt <- names(all_na_alt[!all_na_alt])
  alt_only <- select(alt_only, all_of(all_non_na_alt))

  # Check Numeric Class -----------------------------------------------------
  is_numeric <- apply(alt_only, 2, is.numeric)

  if(!(all(is_numeric))) {
    cli::cli_abort("All alterations in your gene binary must be numeric and only can have values of 0, 1, or NA.
                   Please coerce the following columns to numeric, or pass them to the `other_vars` argument before proceeding: {.field {names(is_numeric[!is_numeric])}}")
  }

  # Check Numeric Class -----------------------------------------------------
  is_numeric <- apply(alt_only, 2, is.numeric)

  if(!(all(is_numeric))) {
    cli::cli_abort("All alterations in your gene binary must be numeric and only can have values of 0, 1, or NA.
                   Please coerce the following columns to numeric, or pass them to the `other_vars` argument before proceeding: {.field {names(is_numeric[!is_numeric])}}")
  }

  # Subset Binary Matrix -----------------------------------------------------
  # get names of genes on the specified panel
  genes_on_panel <- gnomeR::gene_panels %>%
    filter(grepl(panel_id, .data$gene_panel)) %>%
    unnest("genes_in_panel") %>%
    dplyr::pull("genes_in_panel")

  # if not all genes on that panel are included in the genomic binary matrix
  if (length(setdiff(genes_on_panel, names(gene_binary))) != 0){

    cli::cli_alert_info(c("There are {length(setdiff(genes_on_panel, names(gene_binary)))} genes on the {panel_id} panel that are not altered for any patients; columns for those genes are not included in the resulting data frame. ",
    "To see the names of the genes that are not mutated for any patients in the sample, ",
    "To view and compare genes in each panel, see {.code unnest(gnomeR::gene_panels, cols = c('genes_in_panel'))}"))
  }

  # subset to select only those genes
  subset_binary <- select(gene_binary, "sample_id",
                          any_of(other_vars),
                          any_of(genes_on_panel))

  return(subset_binary)
}

