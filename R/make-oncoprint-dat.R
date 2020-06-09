#' dat.oncoPrint
#' Enables creation of a matrix used to generate an OncoPrint heatmap.
#' @param gen.dat A binary matrix or dataframe, with patients as rows and features as columns. Note that the names of the
#' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' @param clin.dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' @return mat : a matrix ready to be plotted using plot.Oncoprint().
#' @export
#' @examples
#' library(gnomeR)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))
#' bin.mut <- binmat(patients = patients,maf = mut,mut.type = "SOMATIC",
#' SNP.only = FALSE,include.silent = FALSE, spe.plat = FALSE)
#' gen.dat <- bin.mut[1:500,
#' names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
#' dat.oncoPrint(gen.dat)
#' ## adding clinical ##
#' clin.patients.dat <-
#' clin.patients[match(abbreviate(rownames(gen.dat),
#' strict = TRUE, minlength = 9),
#' clin.patients$X.Patient.Identifier),] %>%
#' dplyr::rename(DMPID = X.Patient.Identifier,
#'  Smoker = Smoking.History) %>%
#'   dplyr::select(DMPID, Sex,Smoker) %>%
#'   dplyr::filter(!is.na(DMPID)) %>%
#'   dplyr::distinct(DMPID,.keep_all = TRUE)
#' gen.dat <- gen.dat[match(clin.patients.dat$DMPID,
#' abbreviate(rownames(gen.dat),strict = TRUE, minlength = 9)),]
#' clin.patients.dat <- clin.patients.dat %>%
#'   tibble::column_to_rownames('DMPID')
#' rownames(gen.dat) <- rownames(clin.patients.dat)
#' dat.oncoPrint(gen.dat = gen.dat,clin.dat = clin.patients.dat)
#' @import
#' tibble
#' dplyr
#' dtplyr


dat.oncoPrint <- function(gen.dat,clin.dat=NULL){

  # would be best if genetics also had an unknown #
  # if(anyNA(gen.dat))
  #   gen.dat[is.na(gen.dat)] <- 0

  # check if data has non binary cna data and fix it #
  if(length(grep(".cna", colnames(gen.dat))) > 0){
    temp <- gen.dat[,grep(".cna",colnames(gen.dat))]

    temp2 <- do.call("cbind",apply(temp,2,function(x){
      x <- as.numeric(as.character(x))
      yA <- ifelse(x>=0.9,1,0)
      yD <- ifelse(x<=-0.9,1,0)
      out <- as.data.frame(cbind(yA,yD))
      colnames(out) <- c("Amp","Del")
      return(out)
    }))
    colnames(temp2) <- gsub("\\.cna","",colnames(temp2))
    gen.dat <- as.data.frame(cbind(gen.dat[,-grep(".cna",colnames(gen.dat))],
                                   temp2))
  }

  # check that gen.dat is a binary matrix #
  if(anyNA(apply(gen.dat, 2, function(x){anyNA(as.numeric(x[!is.na(x)]))}))){
    warning("All genetic components were not numeric. Those which weren't will be removed")
    gen.dat <- gen.dat[,-which(apply(gen.dat,2,function(x){anyNA(as.numeric(x[!is.na(x)]))}))]
  }

  if(!all(apply(gen.dat,2,function(x){length(unique(x[!is.na(x)])) == 2}))){
    warning("All genetic components were not binary. Those which weren't will be removed")
    gen.dat <- gen.dat[,which(apply(gen.dat,2,function(x){length(unique(x[!is.na(x)])) == 2}))]
  }


  if(!is.null(clin.dat)){
    # seet NA's to UNKNOWN #
    clin.dat <- clin.dat %>%
      tibble::rownames_to_column('sample') %>%
      mutate_all(as.character) %>%
      tibble::column_to_rownames('sample')
    if(anyNA(clin.dat)) clin.dat[is.na(clin.dat)] <- "Unknown"
    # subset data #
    patients <- intersect(rownames(gen.dat),rownames(clin.dat))
    gen.dat <- gen.dat[match(patients,rownames(gen.dat)),]
    clin.dat <- clin.dat %>%
      tibble::rownames_to_column('sample')
    clin.dat <- clin.dat[match(patients,clin.dat$sample),]
    rownames(clin.dat) <- clin.dat$sample
    clin.dat <- clin.dat %>%
      select(-one_of("sample"))
    ## names ##
    genes <- c(colnames(clin.dat),unique(gsub(".Amp|.Del|.fus","",colnames(gen.dat))))
  }
  else genes <- unique(gsub(".Amp|.Del|.fus","",colnames(gen.dat)))

  mat <- as.data.frame(matrix(nrow=length(genes),ncol=nrow(gen.dat)))
  colnames(mat) <- rownames(gen.dat)
  rownames(mat) <- genes

  # mutations #
  mut <- gen.dat[,stats::na.omit(match(genes,colnames(gen.dat)))]

  for(i in 1:nrow(mut)){
    for(j in 1:ncol(mut)){
      if(!is.na(mut[i,j])){
        if(mut[i,j]==1) {
          mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "MUT;"
        }
      }
      else{
        mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "MIS;"
      }
    }
  }

  # del #
  if(length(stats::na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))) > 0){
    del <- gen.dat[,stats::na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))]
    if(is.null(dim(del))) {
      del <- as.data.frame(del)
      rownames(del) <- rownames(gen.dat)
    }
    colnames(del) <- gsub(".Del","",colnames(gen.dat)[stats::na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))])
    for(i in 1:nrow(del)){
      for(j in 1:ncol(del)){
        if(!is.na(del[i,j])){
          if(del[i,j]==1) {
            if(!is.na(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))])){
              mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <-
                paste0(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))],"DEL;")}
            else{mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <- "DEL;"}
          }
        }
      }
    }
  }

  # amp #
  if(length(stats::na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))) > 0){
    amp <- gen.dat[,stats::na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))]
    if(is.null(dim(amp))) {
      amp <- as.data.frame(amp)
      rownames(amp) <- rownames(gen.dat)
    }
    colnames(amp) <- gsub(".Amp","",colnames(gen.dat)[stats::na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))])
    for(i in 1:nrow(amp)){
      for(j in 1:ncol(amp)){
        if(!is.na(amp[i,j])){
          if(amp[i,j]==1) {
            if(!is.na(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))])){
              mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <-
                paste0(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))],"AMP;")}
            else{mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <- "AMP;"}
          }
        }
      }
    }
  }

  # fusions #
  if(length(stats::na.omit(match(paste0(genes,".fus"),colnames(gen.dat)))) > 0){
    fusion <- gen.dat[,stats::na.omit(match(paste0(genes,".fus"),colnames(gen.dat)))]
    if(is.null(dim(fusion))) {
      fusion <- as.data.frame(fusion)
      rownames(fusion) <- rownames(gen.dat)
    }
    colnames(fusion) <- gsub(".fus","",colnames(gen.dat)[stats::na.omit(match(paste0(genes,".fus"),colnames(gen.dat)))])
    for(i in 1:nrow(fusion)){
      for(j in 1:ncol(fusion)){
        if(!is.na(fusion[i,j])){
          if(fusion[i,j]==1) {
            if(!is.na(mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))])){
              mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))] <-
                paste0(mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))],"FUS;")}
            else{mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))] <- "FUS;"}
          }
        }
      }
    }
  }
  ### clin ###
  if(!is.null(clin.dat)){
    temp <- clin.dat
    for(x in colnames(clin.dat)){
      y <- clin.dat[,x]

      # if binary #
      if(length(unique(y[!is.na(y)])) == 2){
        mat[match(x,rownames(mat)),] <- paste0(x,"_",y,";") #ifelse(y==1,"CLIN;",NA)
      }

      else if(length(unique(y[!is.na(y)])) %in% c(3:5)){
        mat[match(x,rownames(mat)),] <- paste0(x,"_",y,";") #abbreviate(x,1)
      }

      else{
        y <- as.numeric(y)
        mat[match(x,rownames(mat)),] <- ifelse(y > stats::median(y,na.rm = T),paste0(x,";"),NA)
      }
    }
  }
  ###
  mat[is.na(mat)] <- "  "
  mat <- as.matrix(mat)

  return(mat)
}
