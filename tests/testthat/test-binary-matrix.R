
# Test Binary Matrix  -----------------------------------------------------------


# Test samples argument ----
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

  # with no fusions in select sample
  gene_binary_no_fusions <-  create_gene_binary(samples=samples,
                                            mutation = gnomeR::mutations,
                                            fusion = gnomeR::sv)

  expect_false(any(str_detect(names(gene_binary_no_fusions), ".fus")))
  expect_equal(nrow(gene_binary_no_fusions), length(samples))
})


test_that("samples selected with no CNA ", {

  samples <- setdiff(gnomeR::mutations$sampleId, gnomeR::cna$sampleId)[1:5]

  # with no fusions in select sample
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


test_that("no samples selected have alterations ", {

  fake_mut <- gnomeR::mutations[1:10, ] %>%
    mutate(sampleId = "a")

  samples <- c("thing", "other_thing")

  expect_error(gene_binary_no_fusions <- create_gene_binary(samples = samples,
                                               mutation = fake_mut,
                                               fusion  = gnomeR::sv), "None*")
})

test_that("some but not all samples selected have alterations ", {

  fake_mut <- gnomeR::mutations[1:10, ] %>%
    mutate(sampleId = "a")

  samples <- c("a", "other_thing")

  expect_no_error(gene_binary <- create_gene_binary(samples = samples, mutation = fake_mut))
  expect_equal(sum(gene_binary[2, 2:4]), 0)
})


# NON UNIQUE SAMPLES in samples ARGUMENT?

# Test data type arguments ------------------------------------------------

# test with and without mut/fusion/cna args passed

# Functions should work with any one of the three passed
test_that("test inputting mut/fusion/cna args can leads to a data.frame output", {

  #Can we obtaine correct result format when either mut/fusion/cna inputted
  expect_error(create_gene_binary())

  expect_true(create_gene_binary( mutation = gnomeR::mutations) %>% is.data.frame())

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


# Test mut_type argument --------------------------------------------------

test_that("test incorrectly specified arg", {

  expect_error(create_gene_binary(mutation = mut2,mut_type = "somatic_only",
                             specify_panel = "no"))
})



# this NA portion was commented out

# test_that("test inclusion of NAs in mut_type ", {
#   mut2 = gnomeR::mutations
#   mut2$mutationStatus[1:10]<-NA
#   mut2$mutationStatus[11:15]<-""
#
#   #example test
#   expect_warning(create_gene_binary(mutation = mut2, specify_panel = "no"), "15 mutations*")
# })



test_that("test inclusion of NAs in mut_type ", {

  mut2 = gnomeR::mutations
  mut2$mutationStatus[1:10]<-NA
  mut2$mutationStatus[11:15]<-""

  # NA included by default (germline_omitted)
  see <- create_gene_binary(mutation = mut2, specify_panel = "no")
  see <- create_gene_binary(mutation = mut2, specify_panel = "no")
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



# Test high_level_cna_only argument --------------------------------------------

test_that("test deletions with -1 and -2 events", {

  test_cna <- tibble::tribble(
    ~hugo_symbol, ~sample_id, ~alteration,
    "TP53",        "samp1",     1,
    "BMPR1A",      "samp2",     -1,
    "CDKN2A",      "samp2",     -1.5,
    "FGFR3",       "samp3",      2,
    "FGFR3",       "samp4",      0,
    "TP53",        "samp2",      -2
  )

  test_mut <- tibble::tribble(
    ~hugo_symbol, ~sample_id,  ~variant_type,  ~mutation_status,   ~variant_classification,
    "TP53",        "samp1",     "SNP",          "Somatic",          "Silent",
    "BMPR1A",      "samp2",     "SNP",          "Somatic",           NA,
    "FGFR3",       "samp3",     "SNP",          "Somatic",           NA,
    "FGFR3",       "samp4",     "SNP",          "Somatic",           NA,
    "TP53",        "samp2",     "SNP",          "Somatic",           NA
  )

  test_fus <- tibble::tribble(
    ~site_1_hugo_symbol, ~site_2_hugo_symbol, ~sample_id,
    "TP53",               "BMPR1A",           "samp1",
    "BMPR1A",              "APC",             "samp2",
    "FGFR3",              "TP53",             "samp3",
    "FGFR3",              "KMT2D",            "samp4",
    "TP53",               NA,                 "samp2"
  )

  samples <- c(test_mut$sample_id,
               test_cna$sample_id,
               test_fus$sample_id) %>% unique()

  proc_all_cna <- create_gene_binary(
    samples = samples,
    mutation = test_mut,
    cna = test_cna,
    fusion = test_fus,
    include_silent = TRUE,
    high_level_cna_only = FALSE)

  proc_hl_cna <- create_gene_binary(
    samples = samples,
    mutation = test_mut,
    cna = test_cna,
    fusion = test_fus,
    include_silent = TRUE,
    high_level_cna_only = TRUE)


  ll_genes <- test_cna %>% filter(alteration == 1 |
                        alteration == -1) %>%
    pull(hugo_symbol)

  diff_gene <- str_remove(
    setdiff(names(proc_all_cna), names(proc_hl_cna)), ".Amp|.Del")

  expect_equal(sort(ll_genes), sort(diff_gene))

})



# Test include_silent argument --------------------------------------------

test_that("test include_silent default when no variant class col", {

  test_cna <- tibble::tribble(
    ~hugo_symbol, ~sample_id, ~alteration,
    "TP53",        "samp1",     1,
    "ALK3",        "samp2",     -1,
    "FGFR3",       "samp3",      2,
    "FGFR3",       "samp4",      0,
    "TP53",        "samp2",      1
  )

  test_mut <- tibble::tribble(
    ~hugo_symbol, ~sample_id,  ~variant_type,  ~mutation_status,
    "TP53",        "samp1",     "SNP",          "Somatic",
    "ALK3",        "samp2",     "SNP",          "Somatic",
    "FGFR3",       "samp3",     "SNP",          "Somatic",
    "FGFR3",       "samp4",     "SNP",          "Somatic",
    "TP53",        "samp2",     "SNP",          "Somatic"
  )

  test_fus <- tibble::tribble(
    ~site_1_hugo_symbol, ~site_2_hugo_symbol, ~sample_id,
    "TP53",               "ALK3",             "samp1",
    "ALK3",               "APC",              "samp2",
    "FGFR3",              "TP53",             "samp3",
    "FGFR3",              "KMT2D",            "samp4",
    "TP53",               NA,                 "samp2"
  )


  samples <- c(test_mut$sample_id,
               test_cna$sample_id,
               test_fus$sample_id) %>% unique()

  expect_error(
    proc <- create_gene_binary(
      samples = samples,
      mutation = test_mut,
      cna = test_cna,
      fusion = test_fus))

  expect_no_error(
    proc <- create_gene_binary(
      mutation = test_mut,
      samples = samples,
      cna = test_cna,
      fusion = test_fus,
      include_silent = TRUE,
      recode_aliases = "no"))

})

test_that("test include_silent silent are removed when variant class col", {

  test_cna <- tibble::tribble(
    ~hugo_symbol, ~sample_id, ~alteration,
    "TP53",        "samp1",     1,
    "BMPR1A",      "samp2",     -1,
    "FGFR3",       "samp3",      2,
    "FGFR3",       "samp4",      0,
    "TP53",        "samp2",      1
  )

  test_mut <- tibble::tribble(
    ~hugo_symbol, ~sample_id,  ~variant_type,  ~mutation_status,   ~variant_classification,
    "TP53",        "samp1",     "SNP",          "Somatic",          "Silent",
    "BMPR1A",      "samp2",     "SNP",          "Somatic",           NA,
    "FGFR3",       "samp3",     "SNP",          "Somatic",           "Other",
    "FGFR3",       "samp4",     "SNP",          "Somatic",           NA,
    "TP53",        "samp2",     "SNP",          "Somatic",           NA
  )

  test_fus <- tibble::tribble(
    ~site_1_hugo_symbol, ~site_2_hugo_symbol, ~sample_id,
    "TP53",               "BMPR1A",             "samp1",
    "BMPR1A",             "APC",              "samp2",
    "FGFR3",              "TP53",             "samp3",
    "FGFR3",              "KMT2D",            "samp4",
    "TP53",               NA,                 "samp2"
  )


  samples <- c(test_mut$sample_id,
               test_cna$sample_id,
               test_fus$sample_id) %>% unique()


  proc_remove_silent1 <- create_gene_binary(
      samples = samples,
      mutation = test_mut,
      cna = test_cna,
      fusion = test_fus)


  proc_remove_silent2 <- create_gene_binary(
    samples = samples,
    mutation = test_mut,
    cna = test_cna,
    include_silent = FALSE,
    fusion = test_fus)

  expect_equal(proc_remove_silent1, proc_remove_silent2)

  proc_keep_silent <-
    create_gene_binary(
      samples = samples,
      mutation = test_mut,
      cna = test_cna,
      include_silent = TRUE,
      fusion = test_fus)

  expect_true(sum(proc_remove_silent1$TP53) <= sum(proc_keep_silent$TP53))

  # make sure FGFGR "Other" type was still included in results
  expect_equal(sum(proc_remove_silent1$FGFR3), sum(proc_keep_silent$FGFR3))

})


# Test snp_only argument --------------------------------------------------

# What happpens  when Variant Type is NA? - Maybe need to add warning to tell user about NAs



# Test the SNP_only argument and recode_aliases argument
test_that("test the snp_only arg", {

  #general tests: input T or F (default is F)
  expect_no_error(gnomeR::create_gene_binary(mutation=gnomeR::mutations, snp_only = T))

  expect_no_warning(gnomeR::create_gene_binary(mutation=gnomeR::mutations,
                                     snp_only = T,
                                     recode_aliases = "no"))
})

#   Test that 0 rows are returned if 0 SNP have been inputted or there is no SNP category in Variant Type been inputted

test_that("test when SNP is NA", {
  mut_snp_zero <- gnomeR::mutations %>%
                 dplyr::filter(.data$variantType != "SNP")

  mut_snp_na <- mut_snp_zero
  mut_snp_na$variantType <- droplevels(factor(mut_snp_na$variantType))

  expect_equal( create_gene_binary(mutation=mut_snp_zero, snp_only = T) %>%
                nrow(),
                0)

  expect_equal( create_gene_binary(mutation=mut_snp_na, snp_only = T) %>%
                  nrow(),
                0)
})
#
#   #What if NA for Variant Type?
#   # note: without Variant Type, the create_gene_binary() still run without error
#   #       snp_only=F will provide full list results and snp_only=T will provide 0 col result
#

test_that("test when SNP is NA", {
  mut_vt_na<- gnomeR::mutations %>%
              dplyr::mutate(variantType=NA)

  expect_true( create_gene_binary(mutation = mut_vt_na, snp_only = F) %>%
               length() >0 )

  expect_equal( create_gene_binary(mutation = mut_vt_na, snp_only = T) %>%
                 nrow(), 0 )
 })



# Other Misc Tests --------------------------------------------------------


