
test_that("Works with binary matrix", {

  bm <- bind_rows(
            "gene10" = c(rep(0, 9), rep(1, 1)),
            "gen50" = c(rep(0, 5), rep(1, 5)),
            "gene20" = c(rep(0, 8), rep(1, 2)),
            "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  # 10% -----
  expect_no_error(sub <- bm %>%
    subset_by_frequency())

  expect_equal(setdiff(names(bm), names(sub)), "gene0")

  # 20% ----
  expect_no_error(sub <- bm %>%
                    subset_by_frequency(t = .2))

  expect_equal(setdiff(names(bm), names(sub)), c("gene10", "gene0"))

  # 50% ----
  expect_no_error(sub <- bm %>%
                    subset_by_frequency(t = .5))

  expect_equal(setdiff(names(bm), names(sub)),
               c("gene10", "gene20", "gene0"))

  # 100% ----
  expect_no_error(sub <- bm %>%
                    subset_by_frequency(t = 1))

  expect_equal(names(sub), "sample_id")

  # 0% ----
  expect_no_error(sub <- bm %>%
                    subset_by_frequency(t = 0))

  expect_equal(sort(names(bm)), sort(names(sub)))

 })


test_that("Error when threshold out of bounds", {

  bm <- bind_rows(
    "gene10" = c(rep(0, 9), rep(1, 1)),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  expect_error(sub <- bm %>%
                    subset_by_frequency(t = -1))

  expect_error(sub <- bm %>%
                 subset_by_frequency(t = 2))

})


test_that("Error when missing sample ID", {

  bm <- bind_rows(
    "gene10" = c(rep(0, 9), rep(1, 1)),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0)))

  expect_error(sub <- bm %>%
                 subset_by_frequency(t = .5))

})


test_that("Error when non numeric cols", {

  bm <- bind_rows(
    "gene10" = c(rep(0, 9), rep(1, 1)),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  bm <- bm %>%
    mutate(across(everything(), ~as.character(.x)))

  expect_error(sub <- bm %>%
                 subset_by_frequency(t = .5))

})

test_that("Counts NAs correctly", {

  bm <- bind_rows(
    "genenow50" = c(0, 0, 0, 0, 0,0,0,1,1,1),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  bm_na <- bind_rows(
    "genenow50" = c(NA, NA, NA, NA, 0,0,0,1,1,1),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  no_na <- bm %>%
    subset_by_frequency(t = .5)

  na <- bm_na %>%
    subset_by_frequency(t = .5)

  expect_equal(setdiff(names(na), names(no_na)), "genenow50")

})

test_that("Correctly removes cols with all NA", {

  bm <- bind_rows(
    "allna" = c(rep(NA, 10)),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  sub <- bm %>%
    subset_by_frequency(t = 0)

  expect_equal(setdiff(names(bm), names(sub)), "allna")

})

