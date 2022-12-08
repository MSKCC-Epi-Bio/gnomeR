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
  pancan = unname(gnomer_colors[c("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                         "GBM",  "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                          "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV",  "PAAD", "PCPG",
                         "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT",
                          "THCA", "THYM", "UCEC", "UCS", "UVM" )]),
  main = unname(gnomer_colors[c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11",
                       "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21","k22",
                       "k23","k24","k25","k26","k27","k28","k29","k30", "k31","k32","k33")]),
  sunset = unname(gnomer_colors[c("k27", "k28", "k29", "k30", "k31", "k32", "k33")])
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
#' @param plot_col Boolean value weather to plot the palette labeled with their hex codes. Defalut is FALSE.
#' @param reverse Boolean indicating whether the palette should be reversed.
#' Default is FALSE.
#' @param ... Additional parameters to pass too `grDevices::colorRampPalette`
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

gnomer_palette <- function(name = "pancan", n, type = c("discrete", "continuous"),
                           plot_col = FALSE, reverse = FALSE,...) {
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

  if (reverse) { pal <- rev(pal)}

  pal_swatch <- switch(type,
         continuous = grDevices::colorRampPalette(pal,...)(n),
         discrete = pal[1:n])

  if(plot_col == FALSE){
    pal_swatch
  }else{
    return(unlist(list(pal_swatch, scales::show_col(pal_swatch))))
  }
}






#' Set gnomeR color palette
#'
#' This function sets the gnomeR color palette as the default palette for all
#' ggplot2 objects. It does so by overriding the following four functions from
#' the ggplot2 package: \code{scale_color_discrete()},
#' \code{scale_fill_discrete()}, \code{scale_color_continuous()}, and
#' \code{scale_fill_continuous()}, and places them in the specified environment.
#' A typical workflow would include this function at the top of a script,
#' and subsequent calls to `ggplot()` will utilize the MSK color palette.
#'
#' @param palette name of palette in gnomer_palettes, supplied in quotes.
#' Options include `"pancan", "main", "sunset"`.
#' Default is `"pancan"`.
#' @param gradient name of gradient palette in `gnomer_palettes`, supplied in quotes.
#' Options include `"pancan", "main", "sunset"`. Default is `"pancan"`.
#' @param reverse if set to `TRUE`, will reverse the order of the color palette
#' @param env environment in which palette will take effect. Default is `rlang::caller_env()`.
#' @author Michael Curry
#' @examples
#' library(ggplot2)
#'
#' set_gnomer_palette()
#'
#' ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point()
#'
#' # setting other MSK palettes
#' set_gnomer_palette(palette = "main", gradient = "sunset")
#'
#' ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point()
#'
#' ggplot(mtcars, aes(wt, mpg, color = cyl)) +
#'   geom_point()
#' @export

set_gnomer_palette <- function(palette = c("pancan", "main", "sunset"),
                               gradient = c("pancan", "main", "sunset"),
                               reverse = FALSE,
                               env = rlang::caller_env()) {
  # choosing the MSK palette
  palette <- gnomer_palettes[[match.arg(palette)]]
  gradient <- gnomer_palettes[[match.arg(gradient)]]

  # reversing palette if requested
  if (reverse) {
    palette <- rev(palette)
    gradient <- rev(gradient)
  }

  # helper functions that set color for continuous ggplot colors
  # both for scales fill and color
  gnomer_fill <- function(...,
                       values = NULL, space = "Lab", na.value = "grey50",
                       guide = "colourbar", aesthetics = "fill") {

    continuous_scale(
      aesthetics = aesthetics,
      scale_name = "gnomercol",
      palette = scales::gradient_n_pal(
        colours = c(gradient[1], gradient[length(gradient)]),
        values = values,
        space = space
      ),
      na.value = na.value,
      guide = guide, ...
    )
  }

  gnomer_colour <- function(...,
                         values = NULL, space = "Lab", na.value = "grey50",
                         guide = "colourbar", aesthetics = "colour") {
    continuous_scale(
      aesthetics = aesthetics,
      scale_name = "gnomercol",
      palette =
        scales::gradient_n_pal(
          colours = c(gradient[1], gradient[length(gradient)]),
          values = values, space = space
        ),
      na.value = na.value,
      guide = guide, ...
    )
  }

  # setting options for ggplot colors
  withr::with_environment(
    env = env,
    code = {
      options("ggplot2.discrete.colour" = unname(palette))
      options("ggplot2.discrete.fill" = unname(palette))
      options("ggplot2.continuous.colour" = gnomer_colour)
      options("ggplot2.continuous.fill" = gnomer_fill)
      options("ggplot2.binned.colour" = unname(palette))
      options("ggplot2.binned.fill" = unname(palette))
    }
  )
}
