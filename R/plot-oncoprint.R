#' plot_oncoPrint
#' Creates the OncoPrint corresponding to the inputted genetic data
#' @param gen.dat A binary matrix or dataframe, with patients as rows and features as columns. Note that the names of the
#' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' @param clin.dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' @param ordered An optional vector of length equal to the number of patients under consideration. Indicates the new order (from left to right)
#' to be plotted.
#' @return p : an oncoprint object
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binmat(patients = patients,maf = mut,
#' fusion = fusion, cna = cna,mut.type = "SOMATIC",
#' SNP.only = FALSE,include.silent = FALSE, spe.plat = FALSE)
#' gen.dat <- bin.mut[,
#' names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
#' plot_oncoPrint(gen.dat)
#' ## adding clinical ##
#' clin.patients.dat <- clin.patients[match(
#' abbreviate(rownames(gen.dat),strict = TRUE, minlength = 9),
#' clin.patients$X.Patient.Identifier),] %>%
#' rename(DMPID = X.Patient.Identifier, Smoker = Smoking.History) %>%
#'   select(DMPID, Sex,Smoker) %>%
#'   filter(!is.na(DMPID)) %>%
#'   distinct(DMPID,.keep_all = TRUE)
#' gen.dat <- gen.dat[match(clin.patients.dat$DMPID,
#' abbreviate(rownames(gen.dat),strict = TRUE,
#'  minlength = 9)),]
#' clin.patients.dat <- clin.patients.dat %>%
#'   tibble::column_to_rownames('DMPID')
#' rownames(gen.dat) <- rownames(clin.patients.dat)
#' plot_oncoPrint(gen.dat = gen.dat,
#' clin.dat = clin.patients.dat)
#' @import
#' ComplexHeatmap
#' tibble

plot_oncoPrint <- function(gen.dat,clin.dat=NULL,ordered=NULL){

  # make data #
  mat <- dat.oncoPrint(gen.dat,clin.dat)
  #############
  if(!is.null(clin.dat)){
    # get all values #
    clin.factors <- unique(unlist(apply(mat,2,unique)))
    clin.factors <- clin.factors[-stats::na.omit(match(c("MUT;","AMP;","DEL;","FUS;","CLIN;","  ",
                                                  "MUT;FUS;","MUT;DEL;","MUT;AMP;"),clin.factors))]
    if(length(clin.factors) == 0){
      col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "CLIN" = "purple")
      alter_fun = list(
        background = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "#CCCCCC", col = NA))
        },
        DEL = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "blue", col = NA))
        },
        AMP = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "red", col = NA))
        },
        MUT = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "#008000", col = NA))
        },
        FUS = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
        },
        CLIN = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "purple", col = NA))
        }
      )

      to.add <- NULL
      added <- NULL
    }

    if(length(clin.factors) > 0){
      clin.factors <- gsub(";","",clin.factors)
      sum.factors <- summary(as.factor(gsub("_.*","",clin.factors)))
      added <- c()
      names.toadd <- c()

      for(k in 1:length(sum.factors)){
        if(as.numeric(sum.factors[k]) == 2){
          added <- c(added, c("black", "#CCCCCC"))
          names.toadd <- c(names.toadd,
                           levels(as.factor(clin.factors[grep(names(sum.factors)[k],clin.factors)])))
        }
        else{
          added <- c(added, grDevices::palette()[c(3,7,2,6,5,4,1)][1:as.numeric(sum.factors[k])])

          names.toadd <- c(names.toadd,
                           levels(as.factor(clin.factors[grep(names(sum.factors)[k],clin.factors)])))
        }
      }
      names(added) <- names.toadd
      col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "CLIN" = "purple",added)
      alter_fun = list(
        background = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "#CCCCCC", col = NA))
        },
        DEL = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "blue", col = NA))
        },
        AMP = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "red", col = NA))
        },
        MUT = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "#008000", col = NA))
        },
        FUS = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
        },
        CLIN = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "purple", col = NA))
        }
      )

      to.add <- lapply(as.character(added),function(val){
        temp <- function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = val, col = NA))
        }
      })
      names(to.add) <- names(added)
      alter_fun <- c(alter_fun,to.add)
    }

    if(!is.null(ordered)) {
      sorted.mat <- mat[,ordered]

      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS",names(added)),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion",names(added))),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP"))))
    }
    else{
      sorted.mat <- mat
      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = NULL,row_order = NULL,
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS",names(added)),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion",names(added))),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP"))))
    }
  }


  if(is.null(clin.dat)){
    col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "CLIN" = "black")
    alter_fun = list(
      background = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "#CCCCCC", col = NA))
      },
      DEL = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "blue", col = NA))
      },
      AMP = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "red", col = NA))
      },
      MUT = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "#008000", col = NA))
      },
      FUS = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
      }
    )

    if(!is.null(ordered)) {
      sorted.mat <- mat[,ordered]

      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS"),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion")),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP"))))
    }
    else{
      sorted.mat <- mat
      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = NULL,row_order = NULL,
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS"),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion")),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP"))))
    }
  }

  return(p)
}

