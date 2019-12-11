#' plot_oncoPrint
#'
#' Creates the OncoPrint corresponding to the inputted genetic data
#'
#' @param gen.dat A binary matrix or dataframe, with patients as rows and features as columns. Note that the names of the
#' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' @param clin.dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' @param ordered An optional vector of length equal to the number of patients under consideration. Indicates the new order (from left to right)
#' to be plotted.
#' @return p : an oncoprint object
#'
#' @export
#'
#' @examples library(gnomeR)
#' mut.only <- create.bin.matrix(maf = mut)
#' all.platforms <- create.bin.matrix(patients = unique(mut$Tumor_Sample_Barcode)[1:500],maf = mut,fusion = fusion,cna = cna)
#' plot_oncoPrint(gen.dat = all.platforms$mut[,c(grep("TP53",colnames(all.platforms$mut)),
#' grep("EGFR",colnames(all.platforms$mut)),
#' grep("KRAS",colnames(all.platforms$mut)),
#' grep("STK11",colnames(all.platforms$mut)),
#' grep("KEAP1",colnames(all.platforms$mut)),
#' grep("ALK",colnames(all.platforms$mut)),
#' grep("CDKN2A",colnames(all.platforms$mut)))])
#' @import ComplexHeatmap


plot_oncoPrint <- function(gen.dat,clin.dat=NULL,ordered=NULL){

  # make data #
  mat <- dat.oncoPrint(gen.dat,clin.dat)
  #############
  col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "CLIN" = "black")
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    DEL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
    },
    AMP = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
    },
    MUT = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
    },
    FUS = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "orange", col = NA))
    },
    CLIN = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
    }
  )


  if(!is.null(ordered)) sorted.mat <- mat[,ordered]
  else sorted.mat <- mat


  if(!is.null(clin.dat)){
  p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                 alter_fun = alter_fun, col = col, column_order = NULL,row_order = NULL,
                 # column_title = "OncoPrint for Meso",
                 barplot_ignore = c("CLIN"),
                 heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","CLIN"),
                                             labels = c("Mutation","Deletion","Amplification","Fusion","ClinicalMarker")))}

  else{
    p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                   alter_fun = alter_fun, col = col, column_order = NULL,row_order = NULL,
                   heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS"),
                                               labels = c("Mutation","Deletion","Amplification","Fusion")))
  }
  return(p)
}




# mat <- dat.oncoPrint(gen.dat = dat.onco[,match(c("BAP1","CDKN2A","CDKN2A.Del","CDKN2B.Del",
#                                                  "NF2","NF2.Del","TP53","TP53.Del","SETD2","TERT"),colnames(dat.onco))],
#                      clin.dat = dat.onco[,match(c("Sex","Stage","Histology","Smoking.Status"),colnames(dat.onco))])
# col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "CLIN" = "black")
# other <- unique(unlist(apply(mat,2,unique)))[-na.omit(match(c("CLIN;","MUT;","DEL;","AMP","  "),unique(unlist(apply(mat,2,unique)))))]
# heatcols <- heat.colors(length(other))
# add.cols <- c()
# for(i in 1:length(heatcols)) {
#   add.cols <- c(add.cols, heatcols[i])
#   names(add.cols)[i] <- gsub(";","",other[i])
# }
# col <- c(col,add.cols)
# alter_fun = list(
#   background = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
#   },
#   DEL = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
#   },
#   AMP = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
#   },
#   MUT = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
#   },
#   FUS = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "orange", col = NA))
#   },
#   CLIN = function(x, y, w, h) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
#   }
# )
# for(i in 1:length(heatcols)){
#   fill.temp <- as.character(add.cols[i])
#   alter_fun[[names(add.cols)[i]]] = function(x, y, w, h,fill = fill.temp) {
#     grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = fill, col = NA)) #as.character(col.temp)
#   }
# }
# sorted.mat <- mat[,mat[3,]=="HBiphasic;"]
# if(!is.null(ordered)) sorted.mat <- mat[,ordered]
# else sorted.mat <- mat
# p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
#                alter_fun = alter_fun, col = col, column_order = NULL,row_order = NULL,
#                barplot_ignore = c("CLIN",names(add.cols)),
#                heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","CLIN",names(add.cols)),
#                                            labels = c("Mutation","Deletion","Amplification","Fusion","ClinicalMarker",names(add.cols))))
