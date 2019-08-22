#' dat.oncoPrint
#'
#' Enables creation of a matrix used to generate an OncoPrint heatmap.
#'
#' @param gen.dat A binary matrix or dataframe, with patients as rows and features as columns. Note that the names of the
#' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' @param clin.dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' @return mat : a matrix ready to be plotted using plot.Oncoprint().
#' @export


dat.oncoPrint <- function(gen.dat,clin.dat=NULL){

  # subset data #
  patients <- intersect(rownames(gen.dat),rownames(clin.dat))
  gen.dat <- gen.dat[match(patients,rownames(gen.dat)),]
  clin.dat <- clin.dat[match(patients,rownames(clin.dat)),]

  ## names ##
  genes <- c(colnames(clin.dat),unique(gsub(".Amp|.Del|.fus","",colnames(gen.dat))))

  mat <- as.data.frame(matrix(nrow=length(genes),ncol=nrow(gen.dat)))
  colnames(mat) <- rownames(gen.dat)
  rownames(mat) <- genes
  # mat[is.na(mat)] <- "  "

  # mutations #
  mut <- gen.dat[,na.omit(match(genes,colnames(gen.dat)))]

  for(i in 1:nrow(mut)){
    for(j in 1:ncol(mut)){
      if(!is.na(mut[i,j])){
        if(mut[i,j]==1) {
          if(rownames(mat)[match(colnames(mut)[j],rownames(mat))] == "WGD") mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "WGD;"
          else{mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "MUT;"}
        }

        if(rownames(mat)[match(colnames(mut)[j],rownames(mat))] == "TMB"){
          if(mut[i,j] < as.numeric(quantile(gen.dat$TMB,0.33)) ){mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "TMB_low;"}
          else if(mut[i,j] > as.numeric(quantile(gen.dat$TMB,0.66)) ){mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "TMB_high;"}
          else{mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "TMB_int;"}
        }

      }
    }
  }

  # del #
  del <- gen.dat[,na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))]
  colnames(del) <- gsub(".Del","",colnames(del))
  for(i in 1:nrow(del)){
    for(j in 1:ncol(del)){
      if(del[i,j]==1) {
        if(!is.na(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))])){
          mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <-
            paste0(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))],"DEL;")}
        else{mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <- "DEL;"}
      }
    }
  }

  # amp #
  amp <- gen.dat[,na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))]
  colnames(amp) <- gsub(".Amp","",colnames(amp))
  for(i in 1:nrow(amp)){
    for(j in 1:ncol(amp)){
      if(amp[i,j]==1) {
        if(!is.na(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))])){
          mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <-
            paste0(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))],"AMP;")}
        else{mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <- "AMP;"}
      }
    }
  }


  ### clin ###

  temp <- clin.dat
  for(x in colnames(clin.dat)){
    y <- clin.dat[,x]

    # if binary #
    if(length(unique(y[!is.na(y)])) == 2){
      mat[match(x,rownames(mat)),] <- ifelse(y==1,"CLIN;",NA)
    }

    else if(length(unique(y[!is.na(y)])) %in% c(3:5)){
      mat[match(x,rownames(mat)),] <- paste0(abbreviate(x,1),y,";")
    }

    else{
      mat[match(x,rownames(mat)),] <- ifelse(y > median(y),"CLIN;",NA)
    }
  }

  ###
  mat[is.na(mat)] <- "  "
  mat <- as.matrix(mat)

  return(mat)
}
