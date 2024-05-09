
test_that("test extracting patient_id from sample_id works", {

  expect_no_error(extract_patient_id(gnomeR::mutations[1:10,]$sampleId))

})

test_that("test error thrown when non IMPACT ID", {

  sample_id =  c("P-0000071-T01-IM3", "XX", "P-0000073-T03-IM5")
  expect_error(extract_patient_id(sample_id))

})

# Test noncolorblindfriendlypairs() function

test_that("test that the non-color blind friendly palettes are properly outputted", {

  expect_no_error(noncolorblindfriendlypairs(pal="main"))

})
