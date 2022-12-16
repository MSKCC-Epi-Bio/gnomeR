# context("check make bin mat")
# source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())
#
# test_that("missing/null variables",{
#   expect_error(binmat(maf=NULL, cna=NULL, fusion=NULL ))
# })
#
# test_that("missing column error",{
#
#
#   mat = mut
#   colnames(mat)[which(colnames(mat)=="Hugo_Symbol")] = "AA"
#   expect_error(binmat(maf=mat), "The MAF file inputted is missing a gene name column. (Hugo_Symbol)", fixed=TRUE)
#
#   mat = mut
#   colnames(mat)[which(colnames(mat)=="Variant_Classification")] = "AA"
#   #expect_error(binmat(maf=mat), "The MAF file inputted is missing a variant classification column. (Variant_Classification)", fixed=TRUE)
#   expect_error(binmat(maf=mat))
#
#   mat = mut
#   colnames(mat)[which(colnames(mat)=="Mutation_Status")] = "AA"
#   expect_warning(binmat(samples=samples, maf =mat ))
#
#   mat = mut
#   colnames(mat)[which(colnames(mat)=="Variant_Type")] = "AA"
#   expect_warning(binmat(samples=samples, maf =mat ))
#
# })
#
# #samples=NULL
#
# test_that("when samples is NULL", {
#   mat = binmat(maf=mut)
#   expect_equal(ncol(mat), 362)
#   expect_equal(nrow(mat), 457)
# })
#
# #specify_panel, also checks samples
# test_that("read in samples with specify_panel", {
#   bin.mut <- binmat(samples = samples,maf = mut,snp_only = F,include_silent = F, specify_panel = T)
#   expect_equal(names(table(colSums(is.na(bin.mut))))[2], "145")
#   expect_equivalent(table(bin.mut[,"TP53"])[2], 124)
# })
#
# #same samples but with different parameters of reading in binmat
#
# test_that("read in samples with snp_only", {
#   #with snp_only=T, some samples have no mutation. this generates a warning.
#   expect_warning(bin.mut <- binmat(samples = samples,maf = mut,snp_only = T,include_silent = F, specify_panel = T))
#
#   #this sample if set.seed is same has no mutations.
#   expect_equal(names(table(t(bin.mut["P-0006282-T02-IM5",]))),"0")
# })
#
#
#

