#' get_tmb
#' \%lifecycle{experimental}
#' Function to calculate the tumor mutation burden of individual patients in a MAF file. Note that this can only be applied to
#' samples sequenced using one of the IMPACT panels. Other samples will be annotated as missing.
#' @param patients a character vector that let's the user specify the patients to be used to create the matrix.
#' Default is NULL is which case all patients in the MAF file will be used.
#' @param maf A MAF file of interest.
#' @param mut.type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param col.names character vector of the necessary columns to be used. By default: col.names = c(Tumor_Sample_Barcode = NULL,
#'  Hugo_Symbol = NULL, Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL).
#' @param ... Further arguments passed belonging to the binmat function.
#' @return TMB A single column dataframe containing the calculated tumor mutation burden of all samples specified adjusted
#'  for the IMPACT panel used for sequencing.
#' @export
#' @examples library(gnomeR)
#' TMB <- get_tmb(maf = mut)
#' @import dplyr
#' @import dtplyr
#' @import stringr

get_tmb <- function(patients = NULL, maf = NULL,
                    mut.type = "SOMATIC",
                    col.names = c(Tumor_Sample_Barcode = NULL, Hugo_Symbol = NULL,
                                  Variant_Classification = NULL, Mutation_Status = NULL, Variant_Type = NULL),
                    ...){
  if(is.null(maf))
    stop("The maf argument must be specified")

  maf <- as_tibble(maf) %>%
    rename(col.names)

  if(is.null(patients))
    patients <- unique(as.character(maf$Tumor_Sample_Barcode))

  if(tolower(mut.type) == "all") Mut.filt = unique(maf$Mutation_Status)
  else Mut.filt = mut.type

  maf <- as_tibble(maf) %>%
    filter(.data$Variant_Classification != "Silent",
           tolower(.data$Mutation_Status) %in% tolower(Mut.filt))

  nb341=sum(ti_341[,3]-ti_341[,2])
  nb410=sum(ti_410[,3]-ti_410[,2])
  nb468=sum(ti_468[,3]-ti_468[,2])

  samples=patients
  v <- strsplit(samples, "-IM|-IH")
  v <- unlist(lapply(1:length(v), function(x)v[[x]][2]))
  idx <- ifelse(v=="3", 1, ifelse(v=="5", 2, ifelse(v=="6", 3, NA)))
  d <- c(nb341, nb410, nb468)

  bb <- NULL
  for(i in 1:length(samples)){
    x=subset(maf, Tumor_Sample_Barcode==samples[i])
    bb[i] <- nrow(x)/d[idx[i]]
  }
  names(bb) <- samples
  bbmb <- bb*1e6

  TMB <- as.data.frame(matrix(nrow = length(samples), ncol = 1))
  rownames(TMB) <- samples
  colnames(TMB) <- "TMB"
  TMB[,1] <- bbmb

  return(TMB)
}

