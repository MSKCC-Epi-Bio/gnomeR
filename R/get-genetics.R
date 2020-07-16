#' Get all genetic data for a study or list of IDs
#'
#' @param sample_ids A character vector of sample ids
#' @param sample_list_id A character vector naming a pre-specified list of
#' samples (e.g. `"mskimpact_Colorectal_Cancer"`)
#' @param genes A list of genes to query. default is all impact genes.
#' @param database A character string of the database to be used. Options are "msk_impact" or "tcga", default is "msk_impact".
#' @param mutations Boolean specifying if mutation data should be fetched. Default is TRUE.
#' @param fusions Boolean specifying if fusion data should be fetched. Default is TRUE.
#' @param cna Boolean specifying if cna data should be fetched. Default is TRUE.
#' @param seg Boolean specifying if the segmentatio data should be fetched. Default is FALSE.
#'
#' @return A dataframe of mutations for each sample ID
#' @export
#'
#' @examples
#' \dontrun{
#' # IMPACT #
#' get_genetics(sample_ids = c("P-0000004-T01-IM3", "P-0000012-T02-IM3"),
#' genes = 207)
#'
#' # TCGA #
#' get_genetics(sample_ids =  c("TCGA-17-Z023-01","TCGA-02-0003-01","TCGA-02-0055-01"),
#'  database = "tcga")
#' }
#' @import
#' cbioportalr

get_genetics <- function(
  sample_ids = NULL,
  sample_list_id = NULL,
  genes = "all_impact",
  database = "msk_impact",
  mutations = TRUE, fusions = TRUE,
  cna = TRUE, seg = FALSE) {

  if(all(c(!mutations, !fusions, !cna, !seg)))
    stop("At least one of the following arguments must be TRUE: mutations, fusions, cna.")
  mut.dat <- NULL
  cna.dat <- NULL
  seg.dat <- NULL


  if(database == "msk_impact"){

  get_cbioportal_db(database)

    if(mutations || fusions){
      mut.dat <- get_mutations(sample_ids,
                              study_id = "mskimpact",
                               genes) %>%
        dplyr::rename(Tumor_Sample_Barcode = "sampleId", Hugo_Symbol = NULL,
               Variant_Classification = "mutationType", Mutation_Status = "mutationStatus",
               Variant_Type = "variantType")
      if(!fusions)
        mut.dat <- mut.dat %>%
          filter(.data$Variant_Classification != "Fusion")
    }

    if(cna)
      cna.dat <- get_cna(sample_ids,
                         study_id = "mskimpact",
                         genes)
    return(list("mut"= mut.dat, "cna" = cna.dat))
  }

  if(database == "tcga"){

    get_cbioportal_db("tcga")

    if(mutations || fusions){
      mut.dat <- get_mutations(sample_ids,
                                    study_id = "all_tcga_studies",
                                    genes) %>%
        dplyr::rename(Tumor_Sample_Barcode = "sampleId", Hugo_Symbol = NULL,
               Variant_Classification = "mutationType", Mutation_Status = "mutationStatus",
               Variant_Type = "variantType")
      if(!fusions)
        mut.dat <- mut.dat %>%
          filter(.data$Variant_Classification != "Fusion")
    }
    if(cna)
      cna.dat <- get_cna(sample_ids,
                              study_id = "all_tcga_studies",
                              genes)
    if(seg)
      seg.dat <- get_segments(sample_ids = sample_ids, study_id = "all_tcga_studies")

    return(list("mut" = mut.dat, "cna" = cna.dat, "seg" = seg.dat))
  }
}
