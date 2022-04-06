# context("check make maf summary")
# source_test_helpers(path = "tests/testthat/helper_dat.R", env = test_env())
#
# test_that("missing column error",{
#
#   mat = mut
#   colnames(mat)[which(colnames(mat)=="Tumor_Sample_Barcode")] = "AA"
#   expect_error(binmat(maf=mat), "The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)", fixed=TRUE)
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
# # test_that("read in 1000 samples with specify.plat", {
# #    mat <- mut %>%
# #     filter(Tumor_Sample_Barcode %in% samples)
# #   plots <- maf.summary(maf = mat)
# #   expect_true(is.ggplot(plots$p.class))
# #   expect_true(is.ggplot(plots$p.type))
# #   expect_true(is.ggplot(plots$p.SNV))
# #   expect_true(is.ggplot(plots$p.patient.variant))
# #   expect_true(is.ggplot(plots$p.variant.bp))
# #   expect_true(is.ggplot(plots$p.genes))
# #   expect_true(is.ggplot(plots$p.variant.dist))
# #   expect_true(is.ggplot(plots$p.variant.dist.bar))
# #   expect_true(is.ggplot(plots$p.SNV.dist))
# #   expect_true(is.ggplot(plots$p.corr))
# #   expect_true(is.ggplot(plots$p.comut))
# #
# #   # expect_is(p$layers[[1]], "proto")
# #   # expect_identical(p$layers[[1]]$geom$objname, "bar")
# #   # expect_identical(p$layers[[1]]$stat$objname, "identity")
# #
# # })
#
#
#
# # test_that("missing column warning but still run",{
# #
# #    mat <- mut %>%
# #     filter(Tumor_Sample_Barcode %in% samples) %>%
# #     select(-one_of("Mutation_Status"))
# #   expect_warning(plots <- maf.summary(maf = mat))
# #   expect_true(is.ggplot(plots$p.class))
# #   expect_true(is.ggplot(plots$p.type))
# #   expect_true(is.ggplot(plots$p.SNV))
# #   expect_true(is.ggplot(plots$p.patient.variant))
# #   expect_true(is.ggplot(plots$p.variant.bp))
# #   expect_true(is.ggplot(plots$p.genes))
# #   expect_true(is.ggplot(plots$p.variant.dist))
# #   expect_true(is.ggplot(plots$p.variant.dist.bar))
# #   expect_true(is.ggplot(plots$p.SNV.dist))
# #   expect_true(is.ggplot(plots$p.corr))
# #   expect_true(is.ggplot(plots$p.comut))
# #
# # })
#
#
# # test_that("Renaming genes",{
# #
# #   mat <- mut %>%
# #     filter(Tumor_Sample_Barcode %in% samples) %>%
# #     mutate(
# #       Hugo_Symbol = as.character(Hugo_Symbol),
# #       Hugo_Symbol = case_when(
# #       Hugo_Symbol == "MLL2" ~ "KMT2D",
# #       Hugo_Symbol == "MLL3" ~ "KMT2C",
# #       TRUE ~ Hugo_Symbol
# #     ))
# #   plots <- maf.summary(maf = mat)
# #   expect_true(is.ggplot(plots$p.class))
# #   expect_true(is.ggplot(plots$p.type))
# #   expect_true(is.ggplot(plots$p.SNV))
# #   expect_true(is.ggplot(plots$p.patient.variant))
# #   expect_true(is.ggplot(plots$p.variant.bp))
# #   expect_true(is.ggplot(plots$p.genes))
# #   expect_true(is.ggplot(plots$p.variant.dist))
# #   expect_true(is.ggplot(plots$p.variant.dist.bar))
# #   expect_true(is.ggplot(plots$p.SNV.dist))
# #   expect_true(is.ggplot(plots$p.corr))
# #   expect_true(is.ggplot(plots$p.comut))
# #
# # })
#
#
# # test_that("Using all mutation statuses",{
# #
# #   mat <- mut %>%
# #     filter(Tumor_Sample_Barcode %in% samples)
# #   plots <- maf.summary(maf = mat,mut.type = "ALL")
# #   expect_true(is.ggplot(plots$p.class))
# #   expect_true(is.ggplot(plots$p.type))
# #   expect_true(is.ggplot(plots$p.SNV))
# #   expect_true(is.ggplot(plots$p.patient.variant))
# #   expect_true(is.ggplot(plots$p.variant.bp))
# #   expect_true(is.ggplot(plots$p.genes))
# #   expect_true(is.ggplot(plots$p.variant.dist))
# #   expect_true(is.ggplot(plots$p.variant.dist.bar))
# #   expect_true(is.ggplot(plots$p.SNV.dist))
# #   expect_true(is.ggplot(plots$p.corr))
# #   expect_true(is.ggplot(plots$p.comut))
# #
# # })
