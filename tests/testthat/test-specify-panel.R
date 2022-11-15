# test specify_impact_panels() ----

test_that("gene binary with muts with impact specified", {
  samples <- gnomeR::cna$sampleId[1:10]
  bin_impact<-  create_gene_binary(samples=samples,
                                           mutation = gnomeR::mutations,
                                           specify_panel = "impact")
  expect_true(ncol(bin_impact) > 0)
  expect_true(nrow(bin_impact) > 1)
})

test_that("gene binary test cna with impact specified", {
  samples <- gnomeR::cna$sampleId[1:10]
  bin_impact<-  create_gene_binary(samples=samples,
                                   cna = gnomeR::cna,
                                   specify_panel = "impact")
  expect_true(ncol(bin_impact) > 0)
  expect_true(nrow(bin_impact) > 1)
})

# test_that("gene binary test fusion with impact specified", {
#   samples <- gnomeR::sv$sampleId[20:30]
#   bin_impact<-  create_gene_binary(samples=samples,
#                                    fusion = gnomeR::sv,
#                                    specify_panel = "impact")
#   expect_true(ncol(bin_impact) > 0)
#   expect_true(nrow(bin_impact) > 1)
# })

test_that("gene binary with all three types of alt and impact only",{
  samples <- gnomeR::mutations$sampleId[gnomeR::mutations$sampleId %in%
                       gnomeR::cna$sampleId]
  samples <-  samples[samples %in% gnomeR::sv$sampleId]
  samples <- samples[1:5]
  bin_impact<-  create_gene_binary(samples=samples,
                                   mutation = gnomeR::mutations,
                                   cna = gnomeR::cna,
                                   fusion = gnomeR::sv)
  expect_true(ncol(bin_impact) > 0)
  expect_true(nrow(bin_impact) > 1)
})

test_that("test 0 impact genes", {
})

