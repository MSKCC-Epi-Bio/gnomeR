
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

test_that("see what happens when no duplicates are in dataset", {


  data <- gnomeR::sv_long[c(1, 20), ]

  expect_no_error(reformat <- reformat_fusion(data))

  expect_equal(nrow(reformat), 2)
  expect_true("TRUE" %in% names(table(reformat[1,] != reformat[2,])))

})

test_that("runs as expected and all fusions remain in dataset", {

  # make all the same person for easy comparison
  # make sure to take off endings of fusion names
  data <- sv_long[1:30, ] %>%
    mutate(sample_id = "TEST")%>%
    mutate(
      #remove leading space in fusion var
      fusion = str_trim(fusion),
      #remove endings to names
      fusion = case_when(
        endsWith(fusion, " fusion") ~ gsub('.{7}$', '', fusion),
        endsWith(fusion,"-intragenic") ~ gsub('.{11}$', "", fusion),
        endsWith(fusion,"-INTRAGENIC") ~ gsub('.{11}$', "", fusion),
        endsWith(fusion,"-INTERGENIC") ~ gsub('.{11}$', "", fusion),
        endsWith(fusion,"-intergenic") ~ gsub('.{11}$', "", fusion),
        endsWith(fusion, " truncation") ~ gsub('.{11}$', "", fusion),
        endsWith(fusion, " rearrangement")  ~ gsub('.{14}$', "", fusion),
        endsWith(fusion, " fusion - Archer") ~ gsub('.{16}$', "", fusion),
        endsWith(fusion, " duplication") ~ gsub('.{11}$', "", fusion),
        endsWith(fusion, " rearrangement") ~ gsub('.{13}$', "", fusion),
        endsWith(fusion, " truncation") ~ gsub('.{10}$', "", fusion),
        endsWith(fusion,  "EZH2(NM_004456) rearrangement exon 5") ~ gsub('.{36}$', "", fusion),
        endsWith(fusion, " PAX5(NM_016734) rearrangement intron 8") ~ gsub('.{38}$', "", fusion),

        TRUE ~ gsub('.{20}$', "", fusion)))

  expect_no_error(reformat <- reformat_fusion(data))

  a <- reformat$event_info
  b <- unique(data$fusion)

  expect_false(length(a) == length(b))

  ################# start here tomorrow #######################

  expect_tru(length(setdiff(b, a)) == length(b) - length(a))

  expect_equal()


})

test <- fusions_sep1 %>% select(c(sample_id, event_info))%>% unique()
save <- setdiff(test, fusions_fus2 %>% select(c(sample_id, event_info)))
