
# Test Binary Matrix Arguments -----------------------------------------------------------

# test patients argument ----
# what happpens when you do specfy it and don't specify it (don't pass anything)?
# what happens when you pass impact samples, non impact samples?


test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})

# test with and without mut/fusion/cna args passed ----
# Functions should work with any one of the three passed
# Does it return results as expected?
# Trying with and without passing patients arg as well- does it return what you'd expect?

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


# test snp_only arg----
# add general tests
# What happpens  when Variant Type is NA? - Maybe need to add warning to tell user about NAs
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})


# test include_silent arg----
# add general tests
# What happens  when Variant_Classification is NA? - Maybe need to add warning to tell user about NAs
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

# test cna_relax arg----
# add general tests
# find an example of data to test this using the API {cbioportalR} that has both 1 and 2 values. Then make it smaller (just a few rows) and test using that
test_that("test", {

  #example test
  expect_equal(TRUE, TRUE)
})


