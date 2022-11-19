
# Test Binary Matrix Arguments -----------------------------------------------------------


# test samples argument ----
# what happens when you pass a vector? What about if you don't specify it (don't pass anything)?
# what happens when you pass impact samples (-IM5/IM6/IM7 etc)?  non impact samples? A mix?


test_that("Check create_gene_binary() subsets based on sample data", {

  mut_valid_sample_ids<-unique(gnomeR::mutations$sampleId)[1:10]

  sub <- create_gene_binary(sample=mut_valid_sample_ids, mutation=gnomeR::mutations)
  expect_message(all <- create_gene_binary(mutation=gnomeR::mutations), "*")

  expect_equal(nrow(sub), length(mut_valid_sample_ids))

  expect_lte(nrow(sub), nrow(all))

})


test_that("Check create_gene_binary() if sample entered in `sample_id` with zero mutations/fusions/cna", {

  mut_valid_sample_ids <- unique(gnomeR::mutations$sampleId)[1:10]

  add_no_mut_sample <- c(mut_valid_sample_ids[1:5], "no_mutations_fake_sample",
                         mut_valid_sample_ids[6:10], "fake")

  gene_binary_no_zero <-  create_gene_binary(sample = mut_valid_sample_ids, mutation = gnomeR::mutations)
  gene_binary_with_zero <-  create_gene_binary(sample = add_no_mut_sample, mutation = gnomeR::mutations)


  expect_equal(gene_binary_with_zero$sample_id, add_no_mut_sample)

  # should be one more obs in data frame with samples arg specified
  expect_equal(nrow(gene_binary_with_zero) -2, length(mut_valid_sample_ids))
})

test_that("samples selected with no fusions ",  {

  samples <- setdiff(gnomeR::mutations$sampleId, gnomeR::sv$sampleId)[1:5]

  # with no fusions in select sample---------
  gene_binary_no_fusions <-  create_gene_binary(samples=samples,
                                            mutation = gnomeR::mutations,
                                            fusion = gnomeR::sv)

  expect_false(any(str_detect(names(gene_binary_no_fusions), ".fus")))
  expect_equal(nrow(gene_binary_no_fusions), length(samples))
})


test_that("samples selected with no CNA ", {

  samples <- setdiff(gnomeR::mutations$sampleId, gnomeR::cna$sampleId)[1:5]

  # with no fusions in select sample---------
  gene_binary_no_cna<-  create_gene_binary(samples=samples,
                                                mutation = gnomeR::mutations,
                                                cna = gnomeR::cna)

  expect_false(any(str_detect(names(gene_binary_no_cna), ".Del")))
  expect_false(any(str_detect(names(gene_binary_no_cna), ".Amp")))

  expect_equal(nrow(gene_binary_no_cna), length(samples))

})

test_that("samples selected with no mutations ", {

  fake_mut <- gnomeR::mutations[1:10, ] %>%
    mutate(sampleId = "a")

  samples <- unique(gnomeR::sv$sampleId[1:5])

  gene_binary_no_fusions <- create_gene_binary(samples = samples,
                                           mutation = fake_mut,
                                           fusion  = gnomeR::sv)

  expect_true(all(str_detect(names(gene_binary_no_fusions)[-1], ".fus")))

  expect_equal(nrow(gene_binary_no_fusions), length(samples))

})


# NON UNIQUE SAMPLES in samples ARGUMENT?

# test with and without mut/fusion/cna args passed ----
# Functions should work with any one of the three passed
# Does it return results as expected?
# Trying with and without passing samples arg as well- does it return what you'd expect?
# what if mut is passed but doesn't have any rows? no columns? no rows or cols?
test_that("test inputting mut/fusion/cna args can leads to a data.frame output", {

  #Can we obtaine correct result format when either mut/fusion/cna inputted
  expect_error(create_gene_binary())

  expect_true( create_gene_binary( mutation = gnomeR::mutations) %>% is.data.frame())

  expect_true( create_gene_binary( fusion = gnomeR::sv) %>% is.data.frame())

  expect_true( create_gene_binary( cna = gnomeR::cna) %>% is.data.frame())

  # What if there is no row/col in the file passing to mutation

  #note: there is no error message if a inputting mut data is 0 rows;
  #      it only return a 0 row and 0 column result
  expect_equal( gnomeR::mutations %>%
                  dplyr::filter(.data$hugoGeneSymbol=="XXXXXXXXX") %>%
                  create_gene_binary(mutation = .) %>%
                  nrow(), 0)

  expect_error( gnomeR::mutations %>%
                dplyr::select(NULL) %>%
                create_gene_binary(mutation = .))


})

# test mut_type argument ----

test_that("test incorrectly specified arg", {

  expect_error(create_gene_binary(mutation = mut2,mut_type = "somatic_only",
                             specify_panel = "no"))
})


test_that("test inclusion of NAs in mut_type ", {
  mut2 = gnomeR::mutations
  mut2$mutationStatus[1:10]<-NA
  mut2$mutationStatus[11:15]<-""

  #example test
  expect_warning(create_gene_binary(mutation = mut2, specify_panel = "no"), "15 mutations*")
})



test_that("test inclusion of NAs in mut_type ", {

  mut2 = gnomeR::mutations
  mut2$mutationStatus[1:10]<-NA
  mut2$mutationStatus[11:15]<-""

  # NA included by default (germline_omitted)
  expect_warning(see <- create_gene_binary(mutation = mut2, specify_panel = "no"))
  check <-see$TP53[which(see$sample_id=="P-0001128-T01-IM3")]
  expect_equal(check, 1)

})

test_that("test inclusion of NAs in mut_type ", {
  mut2 = gnomeR::mutations
  mut2$mutationStatus[1:10]<-NA
  mut2$mutationStatus[11:15]<-""


  # NA included with all
  see = create_gene_binary(mutation = mut2, specify_panel = "no", mut_type = "all")
  expect_equal(see$TP53[which(see$sample_id=="P-0001128-T01-IM3")],1)


  # NA no longer included with somatic_only
  see = create_gene_binary(mutation = mut2, mut_type = "somatic_only", specify_panel = "no")
  expect_equal(see$TP53[which(see$sample_id=="P-0001859-T01-IM3")],0)

  # NA no longer included with germline_only
  see = create_gene_binary(mutation = mut2, mut_type = "germline_only", specify_panel = "no")
  expect_equal(nrow(see), 0)


})


# test specify-panel arg with "impact"- see test_specify-panel.R file



# test snp_only arg----
# add general tests
# What happpens  when Variant Type is NA? - Maybe need to add warning to tell user about NAs
# test_that("test the snp_only arg", {
#
#   #general tests: input T or F (default is F)
#   expect_no_error( create_gene_binary(mutation=gnomeR::mutations, snp_only = T))
#
#   expect_no_warning( create_gene_binary(mutation=gnomeR::mutations,
#                                      snp_only = T,
#                                      recode_aliases = FALSE))
#
#   #TO FIX: this doesn't work with the new example dataset
#   #0 col return if 0 SNP been inputted or there is no SNP category in Variat Type been inputted
#   mut_snp_zero <- gnomeR::mutations %>%
#                  dplyr::filter(.data$variantType != "SNP")
#
#   mut_snp_na <- mut_snp_zero
#   mut_snp_na$variantType <- droplevels(mut_snp_na$variantType)
#
#   expect_equal( create_gene_binary(mutation=mut_snp_zero, snp_only = T) %>%
#                 ncol(),
#                 0)
#
#   expect_equal( create_gene_binary(mutation=mut_snp_na, snp_only = T) %>%
#                   ncol(),
#                 0)
#
#   #What if NA for Variant Type?
#   # note: without Variant Type, the create_gene_binary() still run without error
#   #       snp_only=F will provide full list results and snp_only=T will provide 0 col result
#
#   mut_vt_na<- gnomeR::mutations %>%
#               dplyr::mutate(Variant_Type=NA)
#
#   expect_true( create_gene_binary(mutation = mut_vt_na, snp_only = F) %>%
#                length() >0 )
#
#   expect_equal( create_gene_binary(mutation = mut_vt_na, snp_only = T) %>%
#                  ncol(), 0 )
#
#
# })


# test include_silent arg----
# add general tests
# What happens  when Variant_Classification is NA for some samples in passed data? - Maybe need to add warning to tell user about NAs
# test_that("test include_silent arg", {
#
#   #general tests: input T or F (default is F)
#   expect_error( create_gene_binary(mutation=gnomeR::mutations, include_silent = T), NA)
#
#   expect_warning( create_gene_binary(mutation=gnomeR::mutations,
#                                      include_silent =  T, recode_aliases = FALSE), NA)
#
#
#   #What if NA for variantType?
#   # note: without Variant Type, the create_gene_binary() still run without error
#   #       snp_only=F will provide full list results and snp_only=T will provide 0 col result
#
#   mut_vc_na<- gnomeR::mutations %>%
#     dplyr::mutate(variantType=NA)
#
#   expect_equal( create_gene_binary(mutation = mut_vc_na, include_silent = F) %>% ncol(), 0 )
#
#   expect_true( create_gene_binary(mutation = mut_vc_na, include_silent = T) %>% ncol() > 0 )
#
# })



