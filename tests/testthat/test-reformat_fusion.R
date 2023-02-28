test_that("required columns are included", {
  expect_error(reformat_fusion(gnomeR::sv_long %>% select(-fusion)))
})
