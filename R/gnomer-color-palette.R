#' List of suggested color palettes for when you need a large palette
#'
#' Sometimes you just need a huge palette of fairly distinguishable colors.
#' This is a named vector of Ronglai-approved colors good for things like TCGA PanCan (33 cancer types)
#' or clustering solutions with high K. Run `gnomer_colors` to
#' see the hex codes for the study colors.
#'
#' @export



gnomer_colors <- c(
  #pancan
  "ACC" = "#67000D",
  "BLCA" = "#A50F15",
  "BRCA" = "#EF3B2C",
  "CESC" = "#FC9272",
  "CHOL" = "#FEE0D2",
  "COAD" = "#BCBDDC",
  "DLBC" = "#807DBA",
  "ESCA" = "#54278F",
  "GBM" = "#3F007D",
  "HNSC" = "#08306B",
  "KICH" = "#08519C",
  "KIRC" = "#4292C6",
  "KIRP" = "#9ECAE1",
  "LAML" = "#000000",
  "LGG" = "#525252",
  "LIHC" = "#969696",
  "LUAD" = "#BDBDBD",
  "LUSC" = "#D9D9D9",
  "MESO" = "#80CDC1",
  "OV" = "#35978F",
  "PAAD" = "#01665E",
  "PCPG" = "#006D2C",
  "PRAD" = "#41AB5D",
  "READ" = "#A1D99B",
  "SARC" = "#FFFFCC",
  "SKCM" = "#FED976",
  "STAD" = "#FD8D3C",
  "TGCT" = "#8C510A",
  "THCA" = "#BF812D",
  "THYM" = "#DFC27D",
  "UCEC" = "#FA9FB5",
  "UCS" = "#F768A1",
  "UVM" = "#DD3497",
  # clusters
  "k1" = "#82CC6C",
  "k2" =  "#17A77E",
  "k3" = "#007E7D",
  "k4" =  "#255668",
  "k5" =  "#2A5676B3",
  "k6" =  "#3B809AB3",
  "k7" =  "#67A9B6B3",
  "k8" =  "#99CFD1B3",
  "k9" = "#B7D5E4B3",
  "k10" =  "#8DA3CAB3",
  "k11" =  "#7665A4B3",
  "k12" =  "#6B0077B3",
  "k13" =  "#5B3794B3",
  "k14" =  "#9953A1B3",
  "k15" =  "#C87AADB3",
  "k16" =  "#EBA8BAB3",
  "k17" = "#F8DCD9B3",
  "k18" =  "#E38566B3",
  "k19" =  "#D35D60B3",
  "k20" =  "#B13F63B3",
  "k21" =  "#E24C80B3",
  "k22" =  "#EF7E71B3",
  "k23" =  "#F6A972B3",
  "k24" =  "#FAD18BB3",
  "k25" = "#FDF6B5B3",
  "k26" =  "#DAFF47B3",
  "k27" =  "#EAC800B3",
  "k28" =  "#EB9619B3",
  "k29" =  "#DE655FB3",
  "k30" =  "#C3367EB3",
  "k31" =  "#9C008BB3",
  "k32" =  "#67008CB3",
  "k33" = "#001889B3"
  )



#
# #' Function to extract colors from \code{gnomer_colors} as hex codes
# #'
# #' @param ... Character names of gnomer_colors
# #' @export
#
# gnomer_cols <- function(...) {
#   cols <- c(...)
#
#   if (is.null(cols)) {
#     return(gnomer_colors)
#   }
#
#   unname(gnomer_colors[cols])
# }
#

#' Complete list of gnomeR color palettes
#'
#' @description Creates color palettes based on gnomeR colors, including
#' "main" which is a collection of 33 pale-yet-distinct colors which flow nicely
#' (good for clustering solutions, for example), and
#' "pancan" which is a collection of 33 vibrant and distinct colors, good for visualizing the
#' 33 TCGA PanCan cancer types.
#'
#' @export
#'
#' @examples
#' gnomer_palettes[["pancan"]]

gnomer_palettes <- list(
  `pancan` = gnomer_cols("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                         "GBM",  "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                          "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV",  "PAAD", "PCPG",
                         "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                          "THCA", "THYM", "UCEC", "UCS", "UVM" ),
  `main` = gnomer_cols("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11",
                       "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21","k22",
                       "k23","k24","k25","k26","k27","k28","k29","k30", "k31","k32","k33"),
  `sunset` = gnomer_cols("k27", "k28", "k29", "k30", "k31", "k32", "k33")
)


#' Access the colors in a gnomeR color palette
#'
#' @description gnomeR colors can be accessed and used in plotting
#'
#' @param name Name of desired palette, supplied in quotes. Choices are:
#' "pancan" (default) (best for discrete), "main" (better for discrete), "sunset" (continuous)
#' @param n Number of colors desired. If omitted, uses all colors,
#' or the needed number of colors if less than the total.
#' @param type Either "continuous" or "discrete". Use continuous if you want
#' to automatically interpolate between colours.
#'   @importFrom graphics rgb rect par image text
#'   @importFrom grDevices colorRampPalette
#' @return A vector of colours.
#' @export
#' @keywords colors
#' @examples
#'
#' library(ggplot2)
#'
#' # Print a plot showing the colors in a palette, in order
#'gnomer_palette("pancan")
#'
#' # use a single brand color from a palette
#' # here using the fourth color from the "pancan" palette
#' ggplot(mtcars, aes(hp, mpg)) +
#' geom_point(size = 4, color = gnomer_palette("pancan")[4])
#'
#' # use a discrete color scale - uses fixed colors from the requested palette
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species)) +
#' geom_point(size = 4) +
#' scale_color_manual(values = gnomer_palette("pancan"))
#'
#' # use a continuous color scale - interpolates between colors
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
#' geom_point(size = 4, alpha = .6) +
#' scale_color_gradientn(colors = gnomer_palette("sunset", type = "continuous"))
#'
#' # use a fill color
#' ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
#' geom_histogram(bins = 20, position = "dodge") +
#' scale_fill_manual(values = gnomer_palette("pancan"))
#'

gnomer_palette <- function(name = "pancan", n, type = c("discrete", "continuous")) {
  type <- match.arg(type)

  # since the palettes in msk_palettes are named vectors, ggplot will try to match the names to levels of the variables in the data, which is not what we want. Rather we just want to use them in order, so we need to use unname() here

  pal <- unname(gnomer_palettes[[name]])
  if (is.null(pal)) {
    stop("Palette not found.")
  }

  if (missing(n)) {
    n <- length(pal)
  }

  if (type == "discrete" && n > length(pal)) {
    stop("Number of requested colors greater than what palette can offer")
  }

  out <- switch(type,
                continuous = grDevices::colorRampPalette(pal)(n),
                discrete = pal[1:n]
  )
  structure(out, class = "palette", name = name)
}


#' @importFrom graphics rect par image text
#' @importFrom grDevices rgb
#' @export
print.palette <- function(x, ...) {
  n <- length(x)
  old <- par(mar = c(0.5, 0.5, 0.5, 0.5))
  on.exit(par(old))

  image(1:n, 1, as.matrix(1:n),
        col = x,
        ylab = "", xaxt = "n", yaxt = "n", bty = "n"
  )

  rect(0, 0.9, n + 1, 1.1, col = rgb(1, 1, 1, 0.8), border = NA)
  text((n + 1) / 2, 1, labels = attr(x, "name"), cex = 1, family = "serif")
}


#' Return function to interpolate a gnomeR color palette
#'
#' @param palette Character name of palette in gnomer_palettes.
#' Options include "main", "pancan", "sunset"
#' @param reverse Boolean indicating whether the palette should be reversed.
#' Default to FALSE.
#' @param ... Additional arguments to pass to colorRampPalette()
#' @export

gnomer_pal <- function(palette = "pancan", reverse = FALSE, ...) {
  pal <- gnomer_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  grDevices::colorRampPalette(pal, ...)
}


#' Color scale creator to add gnomeR colors in ggplot
#'
#' @description This color scale generator will interpolate between colors,
#' even when discrete scales are provided.
#' To use exact discrete colors, see examples in \code{gnomer_palette}
#'
#' @param palette Character name of palette in msk_palettes, supplied in quotes.
#' Options include "main" (default), "pancan", "sunset".
#' @param discrete Boolean indicating whether color aesthetic is discrete.
#' Default is TRUE.
#' @param reverse Boolean indicating whether the palette should be reversed.
#' Default is FALSE.
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
#' @export
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # use a discrete color scale
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species)) +
#' geom_point(size = 4) +
#' scale_color_pancan("pancan")
#'
#' # use a continuous color scale
#' ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
#' geom_point(size = 4, alpha = .6) +
#' scale_color_pancan(palette = "sunset", discrete = FALSE)

scale_color_pancan <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- gnomer_pal(palette = palette, reverse = reverse)

  if (discrete) {
    ggplot2::discrete_scale("colour", paste0("gnomer_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_color_gradientn(colours = pal(256), ...)
  }
}


#' Fill scale creator to add gnomeR colors in ggplot
#'
#' @description This fill scale generator will interpolate between the colors
#' in the palette provided
#'
#' @param palette Character name of palette in msk_palettes, supplied in quotes.
#' Options include "main", "pancan", "sunset" (default).
#' @param discrete Boolean indicating whether color aesthetic is discrete or not.
#' Default is TRUE.
#' @param reverse Boolean indicating whether the palette should be reversed.
#' Default is FALSE.
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#' @export
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # use a fill color
#' # alternative use that involves interpolation
#' ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
#' geom_histogram(bins = 20, position = "dodge") +
#' scale_fill_pancan()

scale_fill_pancan <- function(palette = "sunset", discrete = TRUE, reverse = FALSE, ...) {
  pal <- gnomer_pal(palette = palette, reverse = reverse)

  if (discrete) {
    ggplot2::discrete_scale("fill", paste0("gnomer_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_fill_gradientn(colours = pal(256), ...)
  }
}
