#' create.bin.matrix
#'
#' Enables creation of a binary matrix from a maf file with
#' a predifined list of patients (rows are patients and columns are genes)
#'
#' @param patients a character vector that let's the user specify the patients to be used to create the matrix.
#' Default is NULL is which case all patients in the MAF file will be used.
#' @param maf A MAF file.
#' @param mut.type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param SNP.only Boolean to rather the genetics events to be kept only to be SNPs (insertions and deletions will be removed).
#' Default is FALSE.
#' @param include.silent Boolean to keep or remove all silent mutations. TRUE keeps, FALSE removes. Default is FALSE.
#' @return mut : a binary matrix of mutation data
#' @return no.mu.patients : a character vector of patients having no mutations found in the MAF file.
#' @export
#'
#' @import dplyr

create.bin.matrix <- function(patients=NULL, maf, mut.type = "SOMATIC",SNP.only = F,include.silent = F){

  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
  if(length(match("Variant_Classification",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
  if(length(match("Mutation_Status",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a mutation status column. (Mutation_Status)")

  # filter/define patients #
  if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
  else patients <- as.character(unique(maf$Tumor_Sample_Barcode))

  # clean gen dat #
  if(SNP.only) SNP.filt = "SNP" else SNP.filt = unique(maf$Variant_Type)
  if(!include.silent) Variant.filt = "Silent" else Variant.filt = ""
  if(mut.type == "ALL") Mut.filt = unique(maf$Mutation_Status) else Mut.filt = mut.type

  maf <- maf %>% filter(Variant_Classification != Variant.filt,
                        Variant_Type %in% SNP.filt,
                        Mutation_Status %in% Mut.filt)

  #### out frame
  mut <- as.data.frame(matrix(0L,nrow=length(patients),ncol=length(unique(maf$Hugo_Symbol))))
  colnames(mut) <- unique(maf$Hugo_Symbol)
  rownames(mut) <- patients

  for(i in patients){
    genes <- maf$Hugo_Symbol[maf$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0){mut[match(i,rownames(mut)),match(unique(as.character(genes)),colnames(mut))] <- 1}
  }

  missing.mut <- apply(mut,1,function(x){sum(x)==0})
  if(sum(missing.mut) > 0)
    warning("Some patients did not have any mutations found in the MAF file.")

  return(list("mut"=mut,"no.mut.patients"=rownames(mut)[missing.mut]))
}
