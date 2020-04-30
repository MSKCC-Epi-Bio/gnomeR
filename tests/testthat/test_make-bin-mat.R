context("check make bin mat")
source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())

test_that("missing/null variables",{
  expect_error(binmat(maf=NULL, cna=NULL, fusion=NULL ))
})

test_that("missing column error",{

  mat = matrix(rnorm(1,100), ncol=4)
  colnames(mat) = c("Hugo_Symbol", "Variant_Classification", "A","B")
  expect_error(binmat(maf =mat ))

  colnames(mat) = c("A", "Variant_Classification", "Tumor_Sample_Barcode","B")
  expect_error(binmat(maf =mat ))

  colnames(mat) = c("Hugo_Symbol", "A", "Tumor_Sample_Barcode","B")
  expect_error(binmat(maf =mat ))

  mat = mut
  colnames(mat)[which(colnames(mat)=="Mutation_Status")] = "AA"
  patients <- as.character(unique(mat$Tumor_Sample_Barcode))[sample(1:length(unique(mat$Tumor_Sample_Barcode)), 1000, replace=FALSE)]
  expect_warning(binmat(patients=patients, maf =mat ))
})

test_that("read in 1000 patients with spe.plat", {
  bin.mut <- binmat(patients = patients,maf = mut,SNP.only = F,include.silent = F, spe.plat = T, rm.empty = FALSE)
  expect_equal(names(table(colSums(is.na(bin.mut))))[2], "230")
  expect_equivalent(table(bin.mut[,"TP53"])[2], 454 )
})

#same patients but with different parameters of reading in binmat

test_that("read in 1000 patients with SNP.only", {
  #with SNP.only=T, some samples have no mutation. this generates a warning.
  expect_warning(binmat(patients = patients,maf = mut,SNP.only = T,include.silent = F, spe.plat = T, rm.empty = FALSE))
  bin.mut = binmat(patients = patients,maf = mut,SNP.only = T,include.silent = F, spe.plat = T, rm.empty = FALSE)
  #this sample if set.seed is same has no mutations.
  expect_equal(names(table(t(bin.mut["P-0006305-T01-IM5",]))),"0")
})


#test_that("read in 1000 patients with rm.empty", {
  #with SNP.only=T, some samples have no mutation. this generates a warning.
#  bin.mut = binmat(patients = patients,maf = mut,SNP.only = T,include.silent = F, spe.plat = T, rm.empty = TRUE)
  #this sample if set.seed is same has no mutations.
 #expect_equal(nrow(bin.mut),"405")
#})

#and set.plat  - do we need?

#check from 140-145? do we need?





#MAPK3
#PPARG

