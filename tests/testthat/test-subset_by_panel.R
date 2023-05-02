
test_that("Works with binary matrix", {

  bm <- create_gene_binary(mutation = gnomeR::mutations)

  # runs without error
  expect_message(sub <- bm %>%
                    subset_by_panel(panel_id = "IMPACT300"))

  # expect a message since not all genes sequenced on IMPACT300 are mutated
  expect_message(subset_by_panel(bm, panel_id = "IMPACT300"))
})


test_that("Error when panel name incorrectly specified", {
  # panel name not specified
  expect_error(subset_by_panel(bm))

  # panel name is incorrect
  expect_error(subset_by_panel(bm, panel_id = "HACKATHON"))
})


test_that("Error when missing sample ID", {

  expect_error(subset_by_panel(bm %>% select(-sample_id), panel_id = "IMPACT300"))

})


test_that("Error when non numeric cols", {

  bm_sim <- bind_rows(
    "gene10" = c(rep(0, 9), rep(1, 1)),
    "gen50" = c(rep(0, 5), rep(1, 5)),
    "gene20" = c(rep(0, 8), rep(1, 2)),
    "gene0" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.))) %>%
    mutate(across(everything(), ~as.character(.x)))

  expect_error(subset_by_panel(panel = "IMPACT300"))

})

test_that("Correctly removes cols with all NA", {

  bm_sim <- bind_rows(
    "allna" = c(rep(NA, 10)),
    "ALK" = c(rep(0, 5), rep(1, 5))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  sub <- bm_sim %>%
    subset_by_panel(panel_id = "IMPACT300")

  expect_equal(setdiff(names(bm_sim), names(sub)), "allna")

})

test_that("Other columns are retained in `other_vars`", {

  bm_sim <- bind_rows(
    "ALK" = c(rep(0, 5), rep(1, 5)),
    "EGFR" = c(rep(0, 8), rep(1, 2)),
    "BRCA1" = c(rep(0, 10), rep(1, 0)),
    "sex" = rep(c("F", "M"), 5),
    "stage" = rep(c("I", "II"), 5)) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  sub <- bm_sim %>%
    select(-sex, -stage) %>%
    subset_by_frequency(t = 0)

  sub2 <- bm_sim %>%
    subset_by_panel(panel = "IMPACT300", other_vars = c(sex, stage))

  expect_equal(setdiff(names(sub2), names(sub)), c("sex", "stage"))


})

test_that("Pass `other_vars` as strings works", {

  bm_sim <- bind_rows(
    "ALK" = c(rep(0, 5), rep(1, 5)),
    "EGFR" = c(rep(0, 8), rep(1, 2)),
    "BRCA1" = c(rep(0, 10), rep(1, 0)),
    "sex" = rep(c("F", "M"), 5),
    "stage" = rep(c("I", "II"), 5)) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  sub <- bm_sim %>%
    subset_by_panel(panel_id = "IMPACT300", other_vars = c(sex, stage))

  sub2 <- bm_sim %>%
    subset_by_panel(panel_id = "IMPACT300", other_vars = c("sex", "stage"))

  expect_equal(names(sub2), names(sub))
})

test_that("Count of genes not included is correct", {
  bm_sim <- bind_rows(
    "ALK" = c(rep(0, 5), rep(1, 5)),
    "EGFR" = c(rep(0, 8), rep(1, 2)),
    "BRCA1" = c(rep(0, 10), rep(1, 0))) %>%
    mutate(sample_id = as.character(1:nrow(.)))

  expect_message(subset_by_panel(bm_sim, panel_id = "IMPACT300"),
                 "There are 297")
})
