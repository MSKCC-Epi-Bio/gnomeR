context("check make maf summary")

test_that("missing column error",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Variant_Classification", "A","B")
  expect_error(maf.summary(maf = mat ))

})

test_that("missing column warning",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode","B")
  expect_error(maf.summary(maf = mat ))

})
#
# test_that("read in 1000 patients with spe.plat", {
#   set.seed(123)
#   patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
#   bin.mut <- binmat(patients = patients,maf = mut,SNP.only = F,include.silent = F, spe.plat = T, rm.empty = FALSE)
#   expect_equal(names(table(colSums(is.na(bin.mut))))[2], "230")
#   expect_equivalent(table(bin.mut[,"TP53"])[2], 454 )
# })


