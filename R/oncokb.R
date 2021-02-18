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

oncokb <- function(maf = NULL, fusion = NULL, cna = NULL,
                   token = Sys.getenv("ONCOKB_TOKEN"),
                   clin.file = NULL){

  # checks #
  if(is.null(maf) && is.null(fusion) && is.null(cna))
    stop("You must input at least a MAF, fusion or cna file")

  if (identical(token, "")){
  rlang::warn("You must have a valid OncoKB API token, see https://www.oncokb.org.
              You can set a default ONCOKB_TOKEN in `.Renviron` using `usethis::edit_r_environ()`")
    }

  if(is.null(clin.file))
    clin.file = ''

  # check if miniconda is installed + modules #
  # if(!is.character(try(install_miniconda(), silent = T)))
  #   install_miniconda()
  try(install_miniconda(), silent = T)
  # load python files #
  path <- find.package("gnomeR")#.libPaths()[1]#("gnomeR")
  if(is.character(try(source_python(paste0(path,'/AnnotatorCore.py')),silent = TRUE))){
    try(py_install("requests"), silent = T)
    try(py_install("matplotlib"), silent = T)
  }
  source_python(paste0(path,'/AnnotatorCore.py'))
  source_python(paste0(path,'/annotate_maf.py'))
  source_python(paste0(path,'/annotate_fusion.py'))
  source_python(paste0(path,'/annotate_cna.py'))

  # set up #
  maf_oncokb <- NULL
  fusion_oncokb <- NULL
  cna_oncokb <- NULL

  if(!is.null(maf))
    utils::write.table(maf %>%
                  mutate(Protein_position =
                           ifelse(is.na(.data$Protein_position), "", .data$Protein_position)), #%>%
                # select(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, HGVSp)
                file = "temp_maf.txt",quote = F, sep = '\t', row.names = FALSE)
  if(!is.null(fusion))
    write.table(fusion, file = "temp_fusion.txt",quote = F, sep = '\t', row.names = FALSE)
  if(!is.null(cna)){
    if("api" %in% class(cna)){
      temp.cna <- as.data.frame(matrix(0L,
                                       nrow = length(unique(cna$Hugo_Symbol)),
                                       ncol = length(unique(cna$sampleId))+1))
      # rownames(temp.cna) <- unique(cna$SAMPLE_ID)
      temp.cna[,1] <- unique(cna$Hugo_Symbol)
      colnames(temp.cna) <- c("Hugo_Symbol",unique(cna$sampleId))

      for(i in colnames(temp.cna)[-1]){
        temp <- as_tibble(cna) %>%
          filter(.data$sampleId %in% i) %>%
          select(.data$sampleId, .data$Hugo_Symbol, .data$alteration)
        if(nrow(temp)>0){
          temp.cna[match(temp$Hugo_Symbol, temp.cna[,1]),match(i, colnames(temp.cna))] <- temp$alteration
        }
      }
      temp.cna[temp.cna == "Amplification"] <- 2
      temp.cna[temp.cna == "Deletion"] <- -2

      cna <- temp.cna
      temp.cna <- NULL
    }
    write.table(cna, file = "temp_cna.txt",quote = F, sep = '\t', row.names = FALSE)
  }
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
    maf_oncokb <- utils::read.delim('temp_maf_oncoKB.txt')
    file.remove("temp_maf_oncoKB.txt")
  }

  # annotate fusion #
  if(!is.null(fusion)){
    fusionAnnotate(in_fusion = 'temp_fusion.txt', out_fusion = 'temp_fusion_oncoKB.txt',
                   clin_file = clin.file,
                   token = token)
    # remove temp files and load oncokb annotated one #
    file.remove("temp_fusion.txt")
    fusion_oncokb <- read.delim('temp_fusion_oncoKB.txt')
    file.remove("temp_fusion_oncoKB.txt")
  }

  # annotate cna #
  if(!is.null(cna)){

    cnaAnnotate(in_cna = 'temp_cna.txt', out_cna = 'temp_cna_oncoKB.txt',
                clin_file = clin.file,
                token = token)
    # remove temp files and load oncokb annotated one #
    file.remove("temp_cna.txt")
    cna_oncokb <- read.delim('temp_cna_oncoKB.txt')
    file.remove("temp_cna_oncoKB.txt")
  }


  return(list(maf_oncokb = maf_oncokb,
              fusion_oncokb = fusion_oncokb,
              cna_oncokb = cna_oncokb))
}


