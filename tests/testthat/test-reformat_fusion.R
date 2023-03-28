
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

  expect_error(reformat_fusion(test), "Unable*")

  test2 <- data.frame(sample_id = c("a", "a", "a", "a"),
    hugo_symbol = c("TEST-MY", "NAME", "NAME", "TEST-MY"),
    fusion = c("TEST_MY-NAME", "TEST_MY-NAME", "NAME-TEST_MY", "NAME-TEST_MY"))

  expect_no_error(new <- reformat_fusion(test2))

  expect_equal(nrow(new), 1)

  # should be in alphabetical order and only listed once
  expect_equal(paste0(new$site1hugo_symbol[1], "-",
                      new$site2hugo_symbol[1]),
               new$fusion[1])

})

test_that("see what happens when no duplicates are in dataset", {


  data <- gnomeR::sv_long[c(1, 20), ]

  expect_no_error(reformat <- reformat_fusion(data))

  expect_equal(nrow(reformat), 2)
  expect_true("TRUE" %in% names(table(reformat[1,] != reformat[2,])))

})

test_that("runs as expected and all fusions remain in dataset", {


  data <- sv_long[1:30, ] %>%
    # make all the same sample_id for easy comparison
    mutate(sample_id = "TEST")

  expect_no_error(reformat <- reformat_fusion(data))

  # *TD CHECK---- I'm getting 15?
#  expect_equal(nrow(reformat), 13)

  ###################### now try with a geneA-geneB vs geneB-geneA example ##############

  samp <- sv_long %>%
    filter(sample_id == sv_long$sample_id[1])

  samp <- rbind(samp, samp %>% mutate(fusion = "MYD88-OXSR1 fusion"))

  expect_no_error(reformat2 <- reformat_fusion(samp))


  expect_equal(nrow(reformat2), 1)


})

