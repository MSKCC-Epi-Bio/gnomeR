
# Test Plotting Functions-----------------------------------------------------------

#
# # maf_viz --------
#
# test_that("Plot layers match expectations",{
#     samples <- as.character(unique(gnomeR::mut$Tumor_Sample_Barcode))[1:200]
#
#     all_plots <- mutation_viz(mutation = gnomeR::mut %>% filter(Tumor_Sample_Barcode %in% samples))
#
#     expect_equal(names(all_plots), c("varclass",
#                                      "vartype",
#                                      "snvclass",
#                                      "samplevar",
#                                      "topgenes",
#                                      "genecor"))
#     }
# )
#
#
#
# # add.perc --------
# test_that("test add_perc utility", {
#
#   expect_error(ggvarclass(mut) + add.perc(), NA)
# }
#
# )

# gnomeR-color-palette ---------------------------------------------------------
test_that("test incorrrect palette name", {
  expect_error(gnomer_palette("wrong", 4, type = "discrete", plot_col = FALSE,
                              reverse = FALSE), "Palette not found.")
})

test_that("test non-numeric number of colors", {
  expect_error(gnomer_palette("pancan", "wrong", type = "discrete", plot_col = FALSE,
                                        reverse = FALSE), "Specify*")
})

test_that("test missing number of colors", {
  testval <- c("#67000D", "#A50F15", "#EF3B2C", "#FC9272", "#FEE0D2", "#BCBDDC", "#807DBA", "#54278F",
    "#3F007D", "#08306B", "#08519C", "#4292C6", "#9ECAE1", "#000000", "#525252", "#969696",
    "#BDBDBD", "#D9D9D9", "#80CDC1", "#35978F", "#01665E", "#006D2C", "#41AB5D", "#A1D99B",
    "#FFFFCC", "#FED976", "#FD8D3C", "#8C510A", "#BF812D", "#DFC27D", "#FA9FB5", "#F768A1",
    "#DD3497")

  expect_equal(gnomer_palette("pancan", type = "discrete", plot_col = FALSE,
                              reverse = FALSE), testval)
})

test_that("test too many colors requested", {
  expect_error(gnomer_palette("pancan", 35, type = "discrete", plot_col = FALSE,
                              reverse = FALSE), "Number of requested colors greater than what palette can offer")
})

test_that("test reverse option", {
  test_val <- c("#DD3497", "#F768A1", "#FA9FB5", "#DFC27D")

  expect_equal(gnomer_palette("pancan", 4, type = "discrete", plot_col = FALSE,
                              reverse = TRUE), test_val)
})

test_that("test color plots off", {
  test_val <- c("#67000D", "#A50F15", "#EF3B2C", "#FC9272")

  expect_equal(gnomer_palette("pancan", 4, type = "discrete", plot_col = FALSE,
                              reverse = FALSE), test_val)
})

test_that("test color plots on, continuous/discrete", {

  expect_no_error(gnomer_palette("pancan", 4, type = "continuous", plot_col = TRUE,
                              reverse = FALSE))

  expect_no_error(gnomer_palette("pancan", 4, type = "discrete", plot_col = TRUE,
                              reverse = FALSE))
})

# set_gnomer_palette------------------------------------------------------------

test_that("set_gnomer_palette() works", {
  env_palette <- rlang::new_environment()

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         env = env_palette)
      assign("set_gnomer_palette.discrete",
             ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
               geom_point(),
             envir = env_palette)
    }
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         env = env_palette)
      assign("set_gnomer_palette.continuous",
             ggplot(mtcars, aes(wt, mpg, color = cyl)) +
               geom_point(),
             envir = env_palette)
    }
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.discrete",
    get("set_gnomer_palette.discrete", envir = env_palette)
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.continuous",
    get("set_gnomer_palette.continuous", envir = env_palette)
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         reverse = TRUE,
                         env = env_palette)
      assign("set_gnomer_palette.discrete.rev",
             ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
               geom_point(),
             envir = env_palette)
    }
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         reverse = TRUE,
                         env = env_palette)
      assign("set_gnomer_palette.continuous.rev",
             ggplot(mtcars, aes(wt, mpg, color = cyl)) +
               geom_point(),
             envir = env_palette)
    }
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.discrete.rev",
    get("set_gnomer_palette.discrete.rev", envir = env_palette)
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.continuous.rev",
    get("set_gnomer_palette.continuous.rev", envir = env_palette)
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         gradient = "sunset",
                         env = env_palette)
      assign("set_gnomer_palette.discrete.gradient",
             ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
               geom_point(),
             envir = env_palette)
    }
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         gradient = "sunset",
                         env = env_palette)
      assign("set_gnomer_palette.continuous.gradient",
             ggplot(mtcars, aes(wt, mpg, color = cyl)) +
               geom_point(),
             envir = env_palette)
    }
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.discrete.gradient",
    get("set_gnomer_palette.discrete.gradient", envir = env_palette)
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.continuous.gradient",
    get("set_gnomer_palette.continuous.gradient", envir = env_palette)
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         gradient = "sunset",
                         env = env_palette)
      assign("set_gnomer_palette.continuous.gradient.fill",
             ggplot(mtcars, aes(cyl, group = gear, fill = gear)) +
               geom_bar(),
             envir = env_palette)
    }
  )

  vdiffr::expect_doppelganger(
    "set_gnomer_palette.continuous.gradient.fill",
    get("set_gnomer_palette.continuous.gradient.fill", envir = env_palette)
  )

  withr::with_environment(
    env = env_palette,
    code = {
      set_gnomer_palette(palette = "main",
                         env = env_palette)
      reset_gnomer_palette()
      assign("reset_gnomer_palette",
             ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
               geom_point(),
             envir = env_palette)
    }
  )

  vdiffr::expect_doppelganger(
    "reset_gnomer_palette",
    get("reset_gnomer_palette", envir = env_palette)
  )

 }
)


# plotting-functions------------------------------------------------------------
test_that("mutation_viz works", {
  mutations1 <-
    gnomeR::mutations |>
    filter(mutationStatus != "NA")

  mut_plots <- mutation_viz(mutations1)

  varclass_percent <- ggvarclass(mutations1) + add.perc()

  ggcomut_plot <- ggcomut(mutations1)

  vdiffr::expect_doppelganger(
    "mutation_viz_varclass", mut_plots$varclass)

  vdiffr::expect_doppelganger(
    "mutation_viz_varclass_percent", varclass_percent)

  vdiffr::expect_doppelganger(
    "mutation_viz_plots_vartype", mut_plots$vartype)

  vdiffr::expect_doppelganger(
    "mutation_viz_plots_samplevar", mut_plots$samplevar)

  vdiffr::expect_doppelganger(
    "mutation_viz_plots_topgenes", mut_plots$topgenes)

  vdiffr::expect_doppelganger(
    "mutation_viz_plots_genecor", mut_plots$genecor)

  vdiffr::expect_doppelganger(
    "ggcomut_plot", ggcomut_plot)
})

# test_that("varclass plot works", {
#
#   mutation <-
#     mutations %>%
#     rename_columns() %>%
#   mutate(Variant_Classification =
#          stringr::str_replace_all(variant_classification, "_", " ")) %>%
#   mutate(Variant_Classification = variant_classification %>%
#          forcats::fct_infreq() %>%
#          forcats::fct_rev())
#
#   plot1_base <- barplot(height = mutation$Variant_classification,
#                         names = mutation$Variant_classification,
#                         horiz = T)
#   plot1_base
#
#   plot1_ggplot <- mutation_viz(mutations)$varclass
#   vdiffr::expect_doplleganger("varclass", plot1_ggplot)
#
# })

# testthat::test_file("C:\\Users\\toumban\\OneDrive - Memorial Sloan Kettering Cancer Center\\Desktop\\gnomeR\\tests\\testthat\\test-plots.R")

