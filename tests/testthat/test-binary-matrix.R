
# Test Binary Matrix Arguments -----------------------------------------------------------

# General tests ---
test_that("test binary_matrix with mutation runs with no errors", {

  expect_error(binary_matrix(mutation = gnomeR::mut), NA)

  res_mut <- binary_matrix(mutation = gnomeR::mut)
  expect_true(length(res_mut) > 0)
})


test_that("test binary_matrix with cna runs with no errors", {

  expect_error(binary_matrix(cna = gnomeR::cna), NA)

  res <- binary_matrix(cna = gnomeR::cna)
  expect_true(length(res) > 0)
})


test_that("test binary_matrix with fusions runs with no errors", {

  expect_error(binary_matrix(fusion = gnomeR::fusion), NA)

  res_fusion <- binary_matrix(fusion = gnomeR::fusion)
  expect_true(length(res_fusion) > 0)

  length(unique(gnomeR::fusion))
})

test_that("check cna with no alterations are omitted from results", {

  res <- binary_matrix(mutation = gnomeR::mut,
                       cna = gnomeR::cna,
                       fusion = gnomeR::fusion)
  cna_ids <- names(gnomeR::cna)[-1] %>%
    str_replace_all(fixed("."), "-")

  omitted_ids <- setdiff(cna_ids, rownames(res))

  omitted_ids <- omitted_ids %>%
    str_replace_all(fixed("-"), fixed("."))

  check_they_are_zero <- gnomeR::cna %>% select(all_of(omitted_ids)) %>%
    purrr::map_dbl(., ~sum(.x))

  expect_true(sum(check_they_are_zero) == 0)
})

# test samples argument ----
# what happens when you pass a vector? What about if you don't specify it (don't pass anything)?
# what happens when you pass impact samples (-IM5/IM6/IM7 etc)?  non impact samples? A mix?


test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test with and without mut/fusion/cna args passed ----
# Functions should work with any one of the three passed
# Does it return results as expected?
# Trying with and without passing samples arg as well- does it return what you'd expect?
# what if mut is passed but doesn't have any rows? no columns? no rows or cols?
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test mut_type argument ----
# NOTE - SEE EDIT NEEDED IN ISSUE 150 before testing: https://github.com/MSKCC-Epi-Bio/gnomeR/issues/150
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

test_that("test", {
  mut2 = gnomeR::mut
  mut2$Mutation_Status[1:10]<-NA
  mut2$Mutation_Status[11:15]<-""
  #example test
  expect_warning(binary_matrix(mutation = mut2, specify_panel = "no"))
})


test_that("test", {
  mut2 = gnomeR::mut
  mut2$Mutation_Status[1:10]<-NA
  mut2$Mutation_Status[11:15]<-""
  #example test
  see = binary_matrix(mutation = mut2, specify_panel = "no")
  expect_equal(see$TP53[which(rownames(see)=="P-0000062-T01-IM3")],1)
})

test_that("test", {
  mut2 = gnomeR::mut
  mut2$Mutation_Status[1:10]<-NA
  mut2$Mutation_Status[11:15]<-""
  #example test
  see = binary_matrix(mutation = mut2, mut_type = "SOMATIC", specify_panel = "no")
  expect_equal(see$TP53[which(rownames(see)=="P-0000062-T01-IM3")],0)
})


# test snp_only arg----
# add general tests
# What happpens  when Variant Type is NA? - Maybe need to add warning to tell user about NAs
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})


# test include_silent arg----
# add general tests
# What happens  when Variant_Classification is NA for some samples in passed data? - Maybe need to add warning to tell user about NAs
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test cna_binary arg----
# add general tests
# I don't have an example of data that has cna values that aren't just 1 or 2. It would be helpful to
# find an example of data to test this using the API {cbioportalR}. Then make it smaller (just a few rows) and test using that
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test cna_relax arg----
# add general tests
# find an example of data to test this using the API {cbioportalR} that has both 1 and 2 values. Then make it smaller (just a few rows) and test using that
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})



# test rm_empty arg----
# add general tests
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})


