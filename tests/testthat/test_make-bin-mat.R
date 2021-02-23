context("check make bin mat")
source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())

test_that("missing/null variables",{
  expect_error(binmat(maf=NULL, cna=NULL, fusion=NULL ))
})

test_that("missing column error",{

  mat = mut
  colnames(mat)[which(colnames(mat)=="Tumor_Sample_Barcode")] = "AA"
  expect_error(binmat(maf=mat), "The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)", fixed=TRUE)

  mat = mut
  colnames(mat)[which(colnames(mat)=="Hugo_Symbol")] = "AA"
  expect_error(binmat(maf=mat), "The MAF file inputted is missing a gene name column. (Hugo_Symbol)", fixed=TRUE)

  mat = mut
  colnames(mat)[which(colnames(mat)=="Variant_Classification")] = "AA"
  #expect_error(binmat(maf=mat), "The MAF file inputted is missing a variant classification column. (Variant_Classification)", fixed=TRUE)
  expect_error(binmat(maf=mat))

  mat = mut
  colnames(mat)[which(colnames(mat)=="Mutation_Status")] = "AA"
  expect_warning(binmat(patients=patients, maf =mat ))

  mat = mut
  colnames(mat)[which(colnames(mat)=="Variant_Type")] = "AA"
  expect_warning(binmat(patients=patients, maf =mat ))

})

#patients=NULL

test_that("when patients is NULL", {
  mat = binmat(maf=mut)
  expect_equal(ncol(mat), 375)
  expect_equal(nrow(mat), 457)
})

#specify.plat, also checks patients
test_that("read in patients with specify.plat", {
  bin.mut <- binmat(patients = patients,maf = mut,SNP.only = F,include.silent = F, specify.plat = T, rm.empty = FALSE)
  expect_equal(names(table(colSums(is.na(bin.mut))))[2], "145")
  expect_equivalent(table(bin.mut[,"TP53"])[2], 124)
})

#same patients but with different parameters of reading in binmat

test_that("read in patients with SNP.only", {
  #with SNP.only=T, some samples have no mutation. this generates a warning.
  expect_warning(bin.mut <- binmat(patients = patients,maf = mut,SNP.only = T,include.silent = F, specify.plat = T, rm.empty = FALSE))

  #this sample if set.seed is same has no mutations.
  expect_equal(names(table(t(bin.mut["P-0006282-T02-IM5",]))),"0")
})


test_that("read in 300 patients with rm.empty", {
  #with SNP.only=T, some samples have no mutation. this generates a warning.
bin.mut = binmat(patients = patients,maf = mut,SNP.only = T,include.silent = F, specify.plat = T, rm.empty = TRUE)
  #this sample if set.seed is same has no mutations.
  #samples have changed, now there are 457 samples
 expect_equal(ncol(bin.mut),330)
})

#what if we don't have impact samples
test_that("no impact samples",{
  mat = mut
  mat$Tumor_Sample_Barcode = substr(as.character(mat$Tumor_Sample_Barcode),1,12)
  expect_warning(bin.mut<-binmat(maf=mat, specify.plat = TRUE))
  expect_equal(ncol(bin.mut),376)
  #some samples are repeat
  expect_equal(nrow(bin.mut),454)
})



#and set.plat  - do we need?

#check from 140-145? do we need?
#line number 154
#read in fusions





#MAPK3
#PPARG

