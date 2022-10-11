

#thoughts so far:
# still need title case for binary matrix function

cbioportalR::set_cbioportal_db("public")


api_mut <- gnomeR::mut %>%
  mutate(

  )


test_that("test rename_columns runs with no errors", {

  expect_error(rename_columns(gnomeR::mut), NA)

})

# do we want it to display what it changed? currently does not.
test_that("test rename_columns does not have messages", {

  expect_message(rename_columns(gnomeR::mut), NA)

})


test_that("test colnames are renamed properly", {

  expect_snapshot(
    waldo::compare(colnames(gnomeR::mut),
                   colnames(rename_columns(gnomeR::mut)))
  )

})

