#' create.bin.matrix
#'
#' Enables creation of a binary matrix from a maf file with
#' a predifined list of patients (rows are patients and columns are genes)
#'
#' @param patients a character vector that let's the user specify the patients to be used to create the matrix.
#' Default is NULL is which case all patients in the MAF file will be used.
#' @param Fits All fits from the OncoCast run.
#' @return mut : a binary matrix of mutation data
#' @export

create.bin.matrix <- function(patients=NULL, maf){

  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
  if(length(match("Variant_Classification",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")


  if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
  maf <- maf[maf$Variant_Classification != "Silent",]
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
