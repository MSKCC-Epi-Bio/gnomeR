#' OncoKB annotate
#'
#' Enables oncokb annotation of MAF, fusions and CNA files. This is performed using the OncoKB annotator found at https://github.com/oncokb/oncokb-annotator.
#' See details there for file formats.
#'
#' @param maf A maf file to be annotated
#' @param fusion A fusion file to be annotated
#' @param cna A CNA file to be annotated
#' @param token Required token to access OncoKB API, see https://www.oncokb.org/ for details.
#' @param clin.file Optional dataframe containing the cancer types of the samples to be annotated.
#' @return OncoKB annotated files
#' @export
#' @examples
#' \dontrun{
#' library(gnomeR)
#' test <- oncokb(maf = mut[1:100,], token = 'YOUR TOKEN')
#' test$maf_oncokb$oncogenic
#' }
#' @import
#' reticulate
#' dplyr
#' dtplyr

oncokb <- function(maf = NULL, fusion = NULL, cna = NULL, token = '', clin.file = NULL){

  # checks #
  if(is.null(maf) && is.null(fusion) && is.null(cna))
    stop("You must input at least a MAF, fusion or cna file")
  if(token == '')
    stop("You must have a valid token for the OncoKB API, see https://www.oncokb.org/.")
  if(is.null(clin.file))
    clin.file = ''

  # check if miniconda is installed #
  if(!is.character(try(reticulate::install_miniconda(), silent = T)))
    reticulate::install_miniconda()

  # load python files #
  path <- find.package("gnomeR")#.libPaths()[1]#("gnomeR")
  source_python(paste0(path,'/AnnotatorCore.py'))
  source_python(paste0(path,'/annotate_maf.py'))
  source_python(paste0(path,'/annotate_fusion.py'))
  source_python(paste0(path,'/annotate_cna.py'))

  # set up #
  maf_oncokb <- NULL
  fusion_oncokb <- NULL
  cna_oncokb <- NULL

  if(!is.null(maf))
    write.table(maf %>%
                  mutate(Protein_position =
                           ifelse(is.na(Protein_position), "", Protein_position)), #%>%
                # select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, HGVSp)
                file = "temp_maf.txt",quote = F, sep = '\t', row.names = FALSE)
  if(!is.null(fusion))
    write.table(fusion, file = "temp_fusion.txt",quote = F, sep = '\t', row.names = FALSE)
  if(!is.null(cna))
    write.table(cna, file = "temp_cna.txt",quote = F, sep = '\t', row.names = FALSE)
  if(clin.file != ""){
    write.table(clin.file, file = "temp_clin.txt",quote = F, sep = '\t', row.names = FALSE)
    clin.file = "temp_clin.txt"
  }

  # annotate maf #
  if(!is.null(maf)){
    mafAnnotate(in_maf = 'temp_maf.txt', out_maf = 'temp_maf_oncoKB.txt',
                clin_file = clin.file,
                token = token)
    # remove temp files and load oncokb annotated one #
    file.remove("temp_maf.txt")
    maf_oncokb <- read.delim('temp_maf_oncoKB.txt')
    file.remove("temp_maf_oncoKB.txt")
  }

  # annotate fusion #
  if(!is.null(fusion)){
    fusionAnnotate(in_fusion = 'temp_fusion.txt', out_fusion = 'temp_fusion_oncoKB.txt',
                   clin_file = clin.file,
                   token = 'c228f079-e544-4027-a5cd-d2fd3534fd5b')
    # remove temp files and load oncokb annotated one #
    file.remove("temp_fusion.txt")
    fusion_oncokb <- read.delim('temp_fusion_oncoKB.txt')
    file.remove("temp_fusion_oncoKB.txt")
  }

  # annotate cna #
  if(!is.null(cna)){
    cnaAnnotate(in_cna = 'temp_cna.txt', out_cna = 'temp_cna_oncoKB.txt',
                   clin_file = clin.file,
                   token = 'c228f079-e544-4027-a5cd-d2fd3534fd5b')
    # remove temp files and load oncokb annotated one #
    file.remove("temp_cna.txt")
    cna_oncokb <- read.delim('temp_cna_oncoKB.txt')
    file.remove("temp_cna_oncoKB.txt")
  }


  return(list(maf_oncokb = maf_oncokb,
              fusion_oncokb = fusion_oncokb,
              cna_oncokb = cna_oncokb))
}


