context("check make maf summary")

test_that("missing column error (Tumor_Sample_Barcode)",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Variant_Classification", "A","B")
  expect_error(maf.summary(maf = mat ))

})

test_that("missing column error (Hugo_Symbol)",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Tumor_Sample_Barcode", "Variant_Classification", "A","B")
  expect_error(maf.summary(maf = mat ))

})

test_that("missing column warning",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Missing", "Tumor_Sample_Barcode","B")
  expect_error(maf.summary(maf = mat ))
})


test_that("read in 1000 patients with spe.plat", {
  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
  mat <- mut %>%
    filter(Tumor_Sample_Barcode %in% patients)
  plots <- maf.summary(maf = mat)
  expect_true(is.ggplot(plots$p.class))
  expect_true(is.ggplot(plots$p.type))
  expect_true(is.ggplot(plots$p.SNV))
  expect_true(is.ggplot(plots$p.patient.variant))
  expect_true(is.ggplot(plots$p.variant.bp))
  expect_true(is.ggplot(plots$p.genes))
  expect_true(is.ggplot(plots$p.variant.dist))
  expect_true(is.ggplot(plots$p.variant.dist.bar))
  expect_true(is.ggplot(plots$p.SNV.dist))
  expect_true(is.ggplot(plots$p.corr))
  expect_true(is.ggplot(plots$p.comut))

  # expect_is(p$layers[[1]], "proto")
  # expect_identical(p$layers[[1]]$geom$objname, "bar")
  # expect_identical(p$layers[[1]]$stat$objname, "identity")

})



test_that("missing column warning but still run",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
  mat <- mut %>%
    filter(Tumor_Sample_Barcode %in% patients) %>%
    select(-one_of("Mutation_Status"))
  expect_warning(plots <- maf.summary(maf = mat))
  expect_true(is.ggplot(plots$p.class))
  expect_true(is.ggplot(plots$p.type))
  expect_true(is.ggplot(plots$p.SNV))
  expect_true(is.ggplot(plots$p.patient.variant))
  expect_true(is.ggplot(plots$p.variant.bp))
  expect_true(is.ggplot(plots$p.genes))
  expect_true(is.ggplot(plots$p.variant.dist))
  expect_true(is.ggplot(plots$p.variant.dist.bar))
  expect_true(is.ggplot(plots$p.SNV.dist))
  expect_true(is.ggplot(plots$p.corr))
  expect_true(is.ggplot(plots$p.comut))

})


test_that("Renaming genes",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
  mat <- mut %>%
    filter(Tumor_Sample_Barcode %in% patients) %>%
    mutate(
      Hugo_Symbol = as.character(Hugo_Symbol),
      Hugo_Symbol = case_when(
      Hugo_Symbol == "MLL2" ~ "KMT2D",
      Hugo_Symbol == "MLL3" ~ "KMT2C",
      TRUE ~ Hugo_Symbol
    ))
  plots <- maf.summary(maf = mat)
  expect_true(is.ggplot(plots$p.class))
  expect_true(is.ggplot(plots$p.type))
  expect_true(is.ggplot(plots$p.SNV))
  expect_true(is.ggplot(plots$p.patient.variant))
  expect_true(is.ggplot(plots$p.variant.bp))
  expect_true(is.ggplot(plots$p.genes))
  expect_true(is.ggplot(plots$p.variant.dist))
  expect_true(is.ggplot(plots$p.variant.dist.bar))
  expect_true(is.ggplot(plots$p.SNV.dist))
  expect_true(is.ggplot(plots$p.corr))
  expect_true(is.ggplot(plots$p.comut))

})


test_that("Using all mutation statuses",{

  set.seed(123)
  patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
  mat <- mut %>%
    filter(Tumor_Sample_Barcode %in% patients)
  plots <- maf.summary(maf = mat,mut.type = "ALL")
  expect_true(is.ggplot(plots$p.class))
  expect_true(is.ggplot(plots$p.type))
  expect_true(is.ggplot(plots$p.SNV))
  expect_true(is.ggplot(plots$p.patient.variant))
  expect_true(is.ggplot(plots$p.variant.bp))
  expect_true(is.ggplot(plots$p.genes))
  expect_true(is.ggplot(plots$p.variant.dist))
  expect_true(is.ggplot(plots$p.variant.dist.bar))
  expect_true(is.ggplot(plots$p.SNV.dist))
  expect_true(is.ggplot(plots$p.corr))
  expect_true(is.ggplot(plots$p.comut))

})
