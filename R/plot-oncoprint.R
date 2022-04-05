#' plot_oncoPrint
#' Creates the OncoPrint corresponding to the inputted genetic data
#' @param gen.dat A binary matrix or dataframe, with patients as rows and features as columns. Note that the names of the
#' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' @param clin.dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' @param ordered_samples An optional vector of length equal to the number of patients under consideration. Indicates the new order (from left to right)
#' to be plotted.
#' @param rank_genes Boolean argument specifying if the genes should be ranked by event frequency. Default is TRUE.
#' @param background_color Color name to be used to fill the background of the OncoPrint.
#' @return p : an oncoprint object
#' @export
#' @examples library(gnomeR)
#' library(dplyr)
#' library(dtplyr)
#' patients <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#' bin.mut <- binary_matrix(patients = patients, mutation = mut,
#' fusion = fusion, cna = cna,mut.type = "SOMATIC",
#' SNP.only = FALSE,include.silent = FALSE, specify.plat = FALSE)
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
#'   tibble::rownames_to_column("to_rm") %>%
#'   select(-one_of("to_rm")) %>%
#'   tibble::column_to_rownames('DMPID')
#' rownames(gen.dat) <- rownames(clin.patients.dat)
#' plot_oncoPrint(gen.dat = gen.dat,
#' clin.dat = clin.patients.dat)
#' @import
#' ComplexHeatmap
#' tibble

plot_oncoPrint <- function(gen.dat,clin.dat=NULL,ordered_samples=NULL, rank_genes = TRUE,background_color = "#CCCCCC"){


  # if(rank_genes){
  oncoprint_column_order = function(count_matrix = NULL) {
    scoreCol = function(x) {
      score = 0
      for (i in 1:length(x)) {
        if (x[i]) {
          score = score + 2^(length(x) - i * 1/x[i])
        }
      }
      return(score)
    }
    scores = apply(count_matrix[row_order, , drop = FALSE],
                   2, scoreCol)
    order(scores, decreasing = TRUE)
  }
  # }
  # else{
  #   oncoprint_column_order <- function(count_matrix = NULL){
  #     1:ncol(count_matrix)
  #   }
  # }


  # make data #
  mat <- dat.oncoPrint(gen.dat,clin.dat)
  #############
  if(!is.null(clin.dat)){
    # get all values #
    clin.factors <- unique(unlist(apply(mat,2,unique)))
    clin.factors <- clin.factors[-stats::na.omit(match(c("MUT;","AMP;","DEL;","FUS;","  ","MIS;",
                                                         "MUT;FUS;","MUT;DEL;","MUT;AMP;","MUT;DEL;FUS;",
                                                         "DEL;FUS;","MUT;AMP;FUS;",
                                                         "AMP;FUS;"),clin.factors))]
    if(length(clin.factors) == 0){
      col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange","MIS" = "black")
      alter_fun = list(
        background = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = background_color, col = NA))
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
        MIS = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "black", col = NA))
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
          added <- c(added, c("purple", background_color))
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
      col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange","MIS" = "black",added)
      alter_fun = list(
        background = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = background_color, col = NA))
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
        MIS = function(x, y, w, h) {
          grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "black", col = NA))
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

    if(!is.null(ordered_samples)) {
      if(rank_genes)
        row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
      else
        row_order <- 1:nrow(mat)

      sorted.mat <- mat[row_order,ordered_samples]

      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS",names(added)),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion","Missing",names(added))),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
                     alter_fun_is_vectorized = FALSE)
    }
    else{
      if(rank_genes)
        row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
      else
        row_order <- 1:nrow(mat)

      get_type2 = function(x) gsub("^\\s+|\\s+$", "", ComplexHeatmap::default_get_type(x))
      all_type = unique(unlist(lapply(mat, get_type2)))
      all_type = all_type[!is.na(all_type)]
      all_type = all_type[grepl("\\S", all_type)]
      mat_list = lapply(all_type, function(type) {
        m = sapply(mat, function(x) type %in% get_type2(x))
        dim(m) = dim(mat)
        dimnames(m) = dimnames(mat)
        m
      })
      names(mat_list) = all_type[grep("MUT|DEL|AMP|FUS",all_type)]

      arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)),
                  dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
      for (i in seq_along(all_type)) {
        arr[, , i] = mat_list[[i]]
      }
      count_matrix = apply(arr, c(1, 2), sum)

      column_order = oncoprint_column_order(count_matrix = count_matrix)

      sorted.mat <- mat[row_order,column_order]
      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS",names(added)),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion","Missing",names(added))),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
                     alter_fun_is_vectorized = FALSE)
    }
  }


  if(is.null(clin.dat)){
    col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue", "FUS" = "orange", "MIS" = "black")
    alter_fun = list(
      background = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = grid::gpar(fill = background_color, col = NA))
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
      MIS = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = grid::gpar(fill = "black", col = NA))
      }
    )

    if(!is.null(ordered_samples)) {
      if(rank_genes)
        row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
      else
        row_order <- 1:nrow(mat)

      sorted.mat <- mat[row_order,ordered_samples]

      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS"),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion","Missing")),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
                     alter_fun_is_vectorized = FALSE)
    }
    else{
      if(rank_genes)
        row_order <- order(apply(mat, 1, function(x){sum(gsub("MIS;","  ",x) != "  ")}),decreasing = T)
      else
        row_order <- 1:nrow(mat)

      get_type2 = function(x) gsub("^\\s+|\\s+$", "", ComplexHeatmap::default_get_type(x))
      all_type = unique(unlist(lapply(mat, get_type2)))
      all_type = all_type[!is.na(all_type)]
      all_type = all_type[grepl("\\S", all_type)]
      mat_list = lapply(all_type, function(type) {
        m = sapply(mat, function(x) type %in% get_type2(x))
        dim(m) = dim(mat)
        dimnames(m) = dimnames(mat)
        m
      })
      names(mat_list) = all_type[-match("MIS",all_type)]

      arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)),
                  dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
      for (i in seq_along(all_type)) {
        arr[, , i] = mat_list[[i]]
      }
      count_matrix = apply(arr, c(1, 2), sum)

      column_order = oncoprint_column_order(count_matrix = count_matrix)


      sorted.mat <- mat[row_order,column_order]
      p <- oncoPrint(sorted.mat, get_type = function(x) strsplit(x, ";")[[1]],
                     alter_fun = alter_fun, col = col, column_order = 1:ncol(sorted.mat),row_order = 1:nrow(sorted.mat),
                     heatmap_legend_param = list(title = "Alterations", at = c("MUT","DEL","AMP","FUS","MIS"),
                                                 labels = c("Mutation","Deletion","Amplification","Fusion","Missing")),
                     top_annotation = HeatmapAnnotation(
                       column_barplot = anno_oncoprint_barplot(c("MUT","DEL","AMP","FUS"))),
                     alter_fun_is_vectorized = FALSE)
    }
  }

  return(p)
}

