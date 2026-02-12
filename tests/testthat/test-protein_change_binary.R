
# Test Binary Matrix  -----------------------------------------------------------


# Test samples argument ----
# what happens when you pass a vector? What about if you don't specify it (don't pass anything)?
# what happens when you pass impact samples (-IM5/IM6/IM7 etc)?  non impact samples? A mix?


test_that("Check create_protein_change_binary() subsets based on sample data", {

  mut_valid_sample_ids<-unique(gnomeR::mutations$sampleId)[1:10]

  sub <- create_protein_change_binary(samples=mut_valid_sample_ids,
                                      mutation=gnomeR::mutations,
                                      gene=NULL, protein.change=NULL)
  expect_message(all <- create_protein_change_binary(mutation=gnomeR::mutations), "*")

  expect_equal(nrow(sub), length(mut_valid_sample_ids))

  expect_lte(nrow(sub), nrow(all))

})


test_that("Check create_protein_change_binary() if sample entered in `sample_id` not in the mutation file", {

  mut_valid_sample_ids <- unique(gnomeR::mutations$sampleId)[1:10]

  add_no_mut_sample <- c(mut_valid_sample_ids[1:5], "no_mutations_fake_sample",
                         mut_valid_sample_ids[6:10], "fake")

  gene_binary_with_zero <- create_protein_change_binary(samples=add_no_mut_sample,
                                                         mutation = gnomeR::mutations)

  # should be one more obs in data frame with samples arg specified
  expect_equal(nrow(gene_binary_with_zero), length(mut_valid_sample_ids))

  expect_equal(gene_binary_with_zero$sample_id, mut_valid_sample_ids)
})



test_that("no samples selected have mutations", {

  fake_mut <- gnomeR::mutations[1:10, ] %>%
    mutate(sampleId = "a")

  fake_samples <- c("thing", "other_thing")

  expect_error(no_mut <- create_protein_change_binary(samples = fake_samples,
                                                      mutation = fake_mut), "None*")
})


test_that("some but not all samples selected have mutations ", {

  fake_mut <- gnomeR::mutations[1:10, ] %>%
    mutate(sampleId = "a")

  samples <- c("a", "other_thing")

  expect_no_error(protein_change_binary <- create_protein_change_binary(samples = samples,
                                                                        mutation = fake_mut))
  expect_equal(nrow(protein_change_binary), 1)
})



# NON UNIQUE SAMPLES in samples ARGUMENT
test_that("Check when multiple same sample ids are entered", {

  sub_mut = gnomeR::mutations[1:10,]
  samples_dup <- rep(sub_mut$sampleId[1], 2)

  binary_dup <- create_protein_change_binary(samples = samples_dup, mutation = sub_mut)

  expect_equal(binary_dup$sample_id, unique(samples_dup))

})



# Test data type arguments ------------------------------------------------
test_that("Check what happens and what message we get if samples are entered as a dataframe and not a char vect", {

  # get sample IDs that are in both mutations and sv data
  sm_inboth_mf <- merge(x=gnomeR::mutations,y=gnomeR::sv, by="sampleId") %>%
    select(sampleId) %>%
    unique()

  #get fusion data for these IDs
  expect_error(create_protein_change_binary(samples = sm_inboth_mf, fusion =gnomeR::sv))

})



# test with and without mut/fusion/cna args passed

# Functions should work with any one of the three passed
test_that("test inputting mut/fusion/cna args can leads to a data.frame output", {

  #Can we obtain correct result format when either mut/fusion/cna inputted
  expect_error(create_protein_change_binary())

  expect_true(create_protein_change_binary( mutation = gnomeR::mutations) %>% is.data.frame())

})



# Test include_silent argument --------------------------------------------

test_that("test include_silent default when no variant class col", {

  test_mut <- tibble::tribble(
    ~sample_id,          ~hugo_symbol,  ~hgv_sp_short,
    "P-0002296-T01-IM3", "EGFR",	      "L760P",
    "P-0006029-T01-IM5", "EGFR",        "S912A"
  )

  test_binary <- tibble::tribble(
    ~sample_id,          ~EGFR_L760P,  ~EGFR_S912A,
    "P-0002296-T01-IM3",  1,	         0,
    "P-0006029-T01-IM5",  0,           1
  )

  expect_equal(
    create_protein_change_binary(mutation = test_mut,
                                 gene='EGFR',
                                 samples = c("P-0002296-T01-IM3",
                                             "P-0006029-T01-IM5")),
    test_binary
    )
})
