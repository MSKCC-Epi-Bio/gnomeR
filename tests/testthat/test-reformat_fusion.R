
# data checks ---------------------------
test_that("required columns are included & is data.frame", {
  expect_error(reformat_fusion(gnomeR::sv_long %>% select(-fusion)), "The following*")
  expect_error(reformat_fusion(gnomeR::sv_long$hugo_symbol), "`fusion`*")
})

test_that("flags events with two hyphens", {
  test <- data.frame(
    sample_id = c("a", "a", "a", "a"),
    hugo_symbol = c("TEST-MY", "NAME", "NAME", "TEST-MY"),
    fusion = c("TEST-MY-NAME", "TEST-MY-NAME", "NAME-TEST-MY", "NAME-TEST-MY")
  )

  expect_error(reformat_fusion(test), "Some of your hugo_symbols *")

  test2 <- data.frame(sample_id = c("a", "a", "a", "a"),
  hugo_symbol = c("TEST-MY", "NAME", "NAME", "TEST-MY"),
  fusion = c("TEST_MY-NAME", "TEST_MY-NAME", "NAME-TEST_MY", "NAME-TEST_MY"))

  expect_no_error(new <- reformat_fusion(test2))

  expect_equal(nrow(new), 1)

  # should be in alphabetical order and only listed once
  expect_equal(new$event_info[1], "NAME-TEST-MY")
  expect_equal(paste0(new$site1hugo_symbol[1], "-", new$site2hugo_symbol[1]),
               new$event_info[1])

})

test_that("", {

  data <- gnomeR::sv_long[1:30,]
  data2 <- gnomeR::sv_long[c(1, 20), ]

  expect_no_error(reformat <- reformat_fusion(data2))



})

test <- fusions_sep1 %>% select(c(sample_id, event_info))%>% unique()
save <- setdiff(test, fusions_fus2 %>% select(c(sample_id, event_info)))
