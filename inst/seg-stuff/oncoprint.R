#'
#' #' Creates the OncoPrint corresponding to the inputted genetic data
#' #' @param gen_dat A binary matrix or dataframe, with samples as rows and features as columns. Note that the names of the
#' #' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' #' @param clin_dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' #' @param ordered_samples An optional vector of length equal to the number of samples under consideration. Indicates the new order (from left to right)
#' #' to be plotted.
#' #' @param rank_genes Boolean argument specifying if the genes should be ranked by event frequency. Default is TRUE.
#' #' @param background_color Color name to be used to fill the background of the OncoPrint.
#' #' @return an oncoprint object
#' #' @export
#' #' @examples library(gnomeR)
#' #' @import
#' #' ComplexHeatmap
#' #' tibble
#'
#' plot_oncoprint <- function(gen_dat,
#'                            clin_dat=NULL,
#'                            ordered_samples=NULL,
#'                            rank_genes = TRUE,
#'                            background_color = "#CCCCCC"){
#'
#'
#'   # if(rank_genes){
#'   oncoprint_column_order = function(count_matrix = NULL) {
#'     scoreCol = function(x) {
#'       score = 0
#'       for (i in 1:length(x)) {
#'         if (x[i]) {
#'           score = score + 2^(length(x) - i * 1/x[i])
#'         }
#'       }
#'       return(score)
#'     }
#'     scores = apply(count_matrix[row_order, , drop = FALSE],
#'                    2, scoreCol)
#'     order(scores, decreasing = TRUE)
#'   }
#'   # }
#'   # else{
#'   #   oncoprint_column_order <- function(count_matrix = NULL){
#'   #     1:ncol(count_matrix)
#'   #   }
#'   # }
#'
#'
#'   # make data #
#'   mat <- dat_oncoprint(gen_dat,clin_dat)
#'   #############
#'   if(!is.null(clin_dat)){
#'     # get all values #
#'     clin.factors <- unique(unlist(apply(mat,2,unique)))
#'     clin.factors <- clin.factors[-stats::na.omit(match(c("MUT;","AMP;","DEL;","FUS;","  ","MIS;",
#'                                                          "MUT;FUS;","MUT;DEL;","MUT;AMP;","MUT;DEL;FUS;",
#'                                                          "DEL;FUS;","MUT;AMP;FUS;",
#'                                                          "AMP;FUS;"),clin.factors))]
#'     if(length(clin.factors) == 0){
#'       col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange","MIS" = "black")
#'       alter_fun = list(
#'         background = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = background_color, col = NA))
#'         },
#'         DEL = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "blue", col = NA))
#'         },
#'         AMP = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "red", col = NA))
#'         },
#'         MUT = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "#008000", col = NA))
#'         },
#'         FUS = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
#'         },
#'         MIS = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "black", col = NA))
#'         }
#'       )
#'
#'       to.add <- NULL
#'       added <- NULL
#'     }
#'
#'     if(length(clin.factors) > 0){
#'       clin.factors <- gsub(";","",clin.factors)
#'       sum.factors <- summary(as.factor(gsub("_.*","",clin.factors)))
#'       added <- c()
#'       names.toadd <- c()
#'
#'       for(k in 1:length(sum.factors)){
#'         if(as.numeric(sum.factors[k]) == 2){
#'           added <- c(added, c("purple", background_color))
#'           names.toadd <- c(names.toadd,
#'                            levels(as.factor(clin.factors[grep(names(sum.factors)[k],clin.factors)])))
#'         }
#'         else{
#'           added <- c(added, grDevices::palette()[c(3,7,2,6,5,4,1)][1:as.numeric(sum.factors[k])])
#'
#'           names.toadd <- c(names.toadd,
#'                            levels(as.factor(clin.factors[grep(names(sum.factors)[k],clin.factors)])))
#'         }
#'       }
#'       names(added) <- names.toadd
#'       col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange","MIS" = "black",added)
#'       alter_fun = list(
#'         background = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = background_color, col = NA))
#'         },
#'         DEL = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "blue", col = NA))
#'         },
#'         AMP = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "red", col = NA))
#'         },
#'         MUT = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "#008000", col = NA))
#'         },
#'         FUS = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
#'         },
#'         MIS = function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "black", col = NA))
#'         }
#'       )
#'
#'       to.add <- lapply(as.character(added),function(val){
#'         temp <- function(x, y, w, h) {
#'           grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = val, col = NA))
#'         }
#'       })
#'       names(to.add) <- names(added)
#'       alter_fun <- c(alter_fun,to.add)
#'     }
#'
#'     if(!is.null(ordered_samples)) {
#'       if(rank_genes)
#'         row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
#'       else
#'         row_order <- 1:nrow(mat)
#'
#'       sorted.mat <- mat[row_order,ordered_samples]
#'
#'       p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
#'                      alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
#'                      heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS",names(added)),
#'                                                  labels = c("Mutation","Deletion","Amplification","Fusion","Missing",names(added))),
#'                      top_annotation = HeatmapAnnotation(
#'                        column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
#'                      alter_fun_is_vectorized = FALSE)
#'     }
#'     else{
#'       if(rank_genes)
#'         row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
#'       else
#'         row_order <- 1:nrow(mat)
#'
#'       get_type2 = function(x) gsub("^\\s+|\\s+$", "", ComplexHeatmap::default_get_type(x))
#'       all_type = unique(unlist(lapply(mat, get_type2)))
#'       all_type = all_type[!is.na(all_type)]
#'       all_type = all_type[grepl("\\S", all_type)]
#'       mat_list = lapply(all_type, function(type) {
#'         m = sapply(mat, function(x) type %in% get_type2(x))
#'         dim(m) = dim(mat)
#'         dimnames(m) = dimnames(mat)
#'         m
#'       })
#'       names(mat_list) = all_type[grep("MUT|DEL|AMP|FUS",all_type)]
#'
#'       arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)),
#'                   dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
#'       for (i in seq_along(all_type)) {
#'         arr[, , i] = mat_list[[i]]
#'       }
#'       count_matrix = apply(arr, c(1, 2), sum)
#'
#'       column_order = oncoprint_column_order(count_matrix = count_matrix)
#'
#'       sorted.mat <- mat[row_order,column_order]
#'       p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
#'                      alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
#'                      heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS",names(added)),
#'                                                  labels = c("Mutation","Deletion","Amplification","Fusion","Missing",names(added))),
#'                      top_annotation = HeatmapAnnotation(
#'                        column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
#'                      alter_fun_is_vectorized = FALSE)
#'     }
#'   }
#'
#'
#'   if(is.null(clin_dat)){
#'     col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "MIS" = "black")
#'     alter_fun = list(
#'       background = function(x, y, w, h) {
#'         grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = background_color, col = NA))
#'       },
#'       DEL = function(x, y, w, h) {
#'         grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "blue", col = NA))
#'       },
#'       AMP = function(x, y, w, h) {
#'         grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = "red", col = NA))
#'       },
#'       MUT = function(x, y, w, h) {
#'         grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "#008000", col = NA))
#'       },
#'       FUS = function(x, y, w, h) {
#'         grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "orange", col = NA))
#'       },
#'       MIS = function(x, y, w, h) {
#'         grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "black", col = NA))
#'       }
#'     )
#'
#'     if(!is.null(ordered_samples)) {
#'       if(rank_genes)
#'         row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
#'       else
#'         row_order <- 1:nrow(mat)
#'
#'       sorted.mat <- mat[row_order,ordered_samples]
#'
#'       p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
#'                      alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
#'                      heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS"),
#'                                                  labels = c("Mutation","Deletion","Amplification","Fusion","Missing")),
#'                      top_annotation = HeatmapAnnotation(
#'                        column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
#'                      alter_fun_is_vectorized = FALSE)
#'     }
#'     else{
#'       if(rank_genes)
#'         row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
#'       else
#'         row_order <- 1:nrow(mat)
#'
#'       get_type2 = function(x) gsub("^\\s+|\\s+$", "", ComplexHeatmap::default_get_type(x))
#'       all_type = unique(unlist(lapply(mat, get_type2)))
#'       all_type = all_type[!is.na(all_type)]
#'       all_type = all_type[grepl("\\S", all_type)]
#'       mat_list = lapply(all_type, function(type) {
#'         m = sapply(mat, function(x) type %in% get_type2(x))
#'         dim(m) = dim(mat)
#'         dimnames(m) = dimnames(mat)
#'         m
#'       })
#'       names(mat_list) = all_type[-match("MIS",all_type)]
#'
#'       arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)),
#'                   dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
#'       for (i in seq_along(all_type)) {
#'         arr[, , i] = mat_list[[i]]
#'       }
#'       count_matrix = apply(arr, c(1, 2), sum)
#'
#'       column_order = oncoprint_column_order(count_matrix = count_matrix)
#'
#'
#'       sorted.mat <- mat[row_order,column_order]
#'       p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
#'                      alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
#'                      heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS"),
#'                                                  labels = c("Mutation","Deletion","Amplification","Fusion","Missing")),
#'                      top_annotation = HeatmapAnnotation(
#'                        column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
#'                      alter_fun_is_vectorized = FALSE)
#'     }
#'   }
#'
#'   return(p)
#' }
#'
#'
#'
#' #' Enables creation of a matrix used to generate an OncoPrint heatmap.
#' #' @param gen_dat A binary matrix or dataframe, with samples as rows and features as columns. Note that the names of the
#' #' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' #' @param clin_dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' #' @return mat : a matrix ready to be plotted using plot.OncoPrint().
#' #' @export
#' #' @examples
#' #' library(gnomeR)
#' #' library(dplyr)
#' #' samples <- as.character(unique(gnomeR::mutations))
#' #' bin.mut <- create_gene_binary(samples = samples, mutation = gnomeR::mutations,
#' #' mut_type = "somatic_only",
#' #' snp_only = FALSE, include_silent = FALSE, specify_panel = "no")
#' #' gen_dat <- bin.mut[1:500,names(sort(apply(bin.mut,2, sum), decreasing = TRUE))[1:15]]
#' #' dat <- dat_oncoprint(gen_dat)
#' #'
#' #' ## adding clinical ##
#' #' clin.patients.dat <-
#' #' clin.patients[match(abbreviate(rownames(gen_dat),
#' #' strict = TRUE, minlength = 9),
#' #' clin.patients$X.Patient.Identifier),] %>%
#' #' dplyr::rename(DMPID = X.Patient.Identifier,
#' #'  Smoker = Smoking.History) %>%
#' #'   dplyr::select(DMPID, Sex,Smoker) %>%
#' #'   dplyr::filter(!is.na(DMPID)) %>%
#' #'   dplyr::distinct(DMPID,.keep_all = TRUE)
#' #' gen_dat <- gen_dat[match(clin.patients.dat$DMPID,
#' #' abbreviate(rownames(gen_dat),strict = TRUE, minlength = 9)),]
#' #' clin.patients.dat <- clin.patients.dat %>%
#' #' dplyr::mutate(DMPID = as.character(DMPID)) %>%
#' #'   tibble::rownames_to_column("to_rm") %>%
#' #'   select(-one_of("to_rm")) %>%
#' #'   tibble::column_to_rownames('DMPID')
#' #' rownames(gen_dat) <- rownames(clin.patients.dat)
#' #' dat <- dat_oncoprint(gen_dat = gen_dat, clin_dat = clin.patients.dat)
#' #' @import
#' #' tibble
#' #' dplyr
#'
#'
#' dat_oncoprint <- function(gen_dat, clin_dat=NULL){
#'
#'   # would be best if genetics also had an unknown #
#'   # if(anyNA(gen_dat))
#'   #   gen_dat[is.na(gen_dat)] <- 0
#'
#'   # check if data has non binary cna data and fix it #
#'   if(length(grep(".cna", colnames(gen_dat))) > 0){
#'     temp <- gen_dat[,grep(".cna",colnames(gen_dat))]
#'
#'     temp2 <- do.call("cbind",apply(temp,2,function(x){
#'       x <- as.numeric(as.character(x))
#'       yA <- ifelse(x>=0.9,1,0)
#'       yD <- ifelse(x<=-0.9,1,0)
#'       out <- as.data.frame(cbind(yA,yD))
#'       colnames(out) <- c("Amp","Del")
#'       return(out)
#'     }))
#'     colnames(temp2) <- gsub("\\.cna","",colnames(temp2))
#'     gen_dat <- as.data.frame(cbind(gen_dat[,-grep(".cna",colnames(gen_dat))],
#'                                    temp2))
#'   }
#'
#'   # check that gen_dat is a binary matrix #
#'   if(anyNA(apply(gen_dat, 2, function(x){anyNA(as.numeric(x[!is.na(x)]))}))){
#'     warning("All genetic components were not numeric. Those which weren't will be removed")
#'     gen_dat <- gen_dat[,-which(apply(gen_dat,2,function(x){anyNA(as.numeric(x[!is.na(x)]))}))]
#'   }
#'
#'   if(!all(apply(gen_dat,2,function(x){length(unique(x[!is.na(x)])) <= 2}))){
#'     warning("All genetic components were not binary. Those which weren't will be removed")
#'     gen_dat <- gen_dat[,-which(apply(gen_dat,2,function(x){length(unique(x[!is.na(x)])) > 2}))]
#'   }
#'
#'
#'   if(!is.null(clin_dat)){
#'     # seet NA's to UNKNOWN #
#'     clin_dat <- clin_dat %>%
#'       tibble::rownames_to_column('sample') %>%
#'       mutate_all(as.character) %>%
#'       tibble::column_to_rownames('sample')
#'     if(anyNA(clin_dat)) clin_dat[is.na(clin_dat)] <- "Unknown"
#'     # subset data #
#'     patients <- intersect(rownames(gen_dat),rownames(clin_dat))
#'     gen_dat <- gen_dat[match(patients,rownames(gen_dat)),]
#'     clin_dat <- clin_dat %>%
#'       tibble::rownames_to_column('sample')
#'     clin_dat <- clin_dat[match(patients,clin_dat$sample),]
#'     rownames(clin_dat) <- clin_dat$sample
#'     clin_dat <- clin_dat %>%
#'       select(-one_of("sample"))
#'     ## names ##
#'     genes <- c(colnames(clin_dat),unique(gsub(".Amp|.Del|.fus","",colnames(gen_dat))))
#'   }
#'   else genes <- unique(gsub(".Amp|.Del|.fus","",colnames(gen_dat)))
#'
#'   mat <- as.data.frame(matrix(nrow=length(genes),ncol=nrow(gen_dat)))
#'   colnames(mat) <- rownames(gen_dat)
#'   rownames(mat) <- genes
#'
#'   # mutations #
#'   mut <- gen_dat[,stats::na.omit(match(genes,colnames(gen_dat)))]
#'
#'   for(i in 1:nrow(mut)){
#'     for(j in 1:ncol(mut)){
#'       if(!is.na(mut[i,j])){
#'         if(mut[i,j]==1) {
#'           mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "MUT;"
#'         }
#'       }
#'       else{
#'         mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "MIS;"
#'       }
#'     }
#'   }
#'
#'   # del #
#'   if(length(stats::na.omit(match(paste0(genes,".Del"),colnames(gen_dat)))) > 0){
#'     del <- gen_dat[,stats::na.omit(match(paste0(genes,".Del"),colnames(gen_dat)))]
#'     if(is.null(dim(del))) {
#'       del <- as.data.frame(del)
#'       rownames(del) <- rownames(gen_dat)
#'     }
#'     colnames(del) <- gsub(".Del","",colnames(gen_dat)[stats::na.omit(match(paste0(genes,".Del"),colnames(gen_dat)))])
#'     for(i in 1:nrow(del)){
#'       for(j in 1:ncol(del)){
#'         if(!is.na(del[i,j])){
#'           if(del[i,j]==1) {
#'             if(!is.na(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))])){
#'               mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <-
#'                 paste0(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))],"DEL;")}
#'             else{mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <- "DEL;"}
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   # amp #
#'   if(length(stats::na.omit(match(paste0(genes,".Amp"),colnames(gen_dat)))) > 0){
#'     amp <- gen_dat[,stats::na.omit(match(paste0(genes,".Amp"),colnames(gen_dat)))]
#'     if(is.null(dim(amp))) {
#'       amp <- as.data.frame(amp)
#'       rownames(amp) <- rownames(gen_dat)
#'     }
#'     colnames(amp) <- gsub(".Amp","",colnames(gen_dat)[stats::na.omit(match(paste0(genes,".Amp"),colnames(gen_dat)))])
#'     for(i in 1:nrow(amp)){
#'       for(j in 1:ncol(amp)){
#'         if(!is.na(amp[i,j])){
#'           if(amp[i,j]==1) {
#'             if(!is.na(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))])){
#'               mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <-
#'                 paste0(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))],"AMP;")}
#'             else{mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <- "AMP;"}
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   # fusions #
#'   if(length(stats::na.omit(match(paste0(genes,".fus"),colnames(gen_dat)))) > 0){
#'     fusion <- gen_dat[,stats::na.omit(match(paste0(genes,".fus"),colnames(gen_dat)))]
#'     if(is.null(dim(fusion))) {
#'       fusion <- as.data.frame(fusion)
#'       rownames(fusion) <- rownames(gen_dat)
#'     }
#'     colnames(fusion) <- gsub(".fus","",colnames(gen_dat)[stats::na.omit(match(paste0(genes,".fus"),colnames(gen_dat)))])
#'     for(i in 1:nrow(fusion)){
#'       for(j in 1:ncol(fusion)){
#'         if(!is.na(fusion[i,j])){
#'           if(fusion[i,j]==1) {
#'             if(!is.na(mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))])){
#'               mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))] <-
#'                 paste0(mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))],"FUS;")}
#'             else{mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))] <- "FUS;"}
#'           }
#'         }
#'       }
#'     }
#'   }
#'   ### clin ###
#'   if(!is.null(clin_dat)){
#'     temp <- clin_dat
#'     for(x in colnames(clin_dat)){
#'       y <- clin_dat[,x]
#'
#'       # if binary #
#'       if(length(unique(y[!is.na(y)])) == 2){
#'         mat[match(x,rownames(mat)),] <- paste0(x,"_",y,";") #ifelse(y==1,"CLIN;",NA)
#'       }
#'
#'       else if(length(unique(y[!is.na(y)])) %in% c(3:5)){
#'         mat[match(x,rownames(mat)),] <- paste0(x,"_",y,";") #abbreviate(x,1)
#'       }
#'
#'       else{
#'         y <- as.numeric(y)
#'         mat[match(x,rownames(mat)),] <- ifelse(y > stats::median(y,na.rm = T),paste0(x,";"),NA)
#'       }
#'     }
#'   }
#'   ###
#'   mat[is.na(mat)] <- "  "
#'   mat <- as.matrix(mat)
#'
#'   return(mat)
#' }
#'
#'
