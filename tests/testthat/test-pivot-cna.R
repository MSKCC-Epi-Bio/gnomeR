
# Pivot Wider -----------------------------------------------------------------

test_that("pivot CNA- check amplification numbers", {

  cna_long <- data.frame(
      sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                   "P-0005436-T01-IM3",
                   "P-0001276-T01-IM3","P-0003333-T01-IM3"),
      Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                      "HIST1H3B","KDR"),
      alteration = c("AMPLIFICATION","AMPLIFICATION",
                     "AMPLIFICATION","AMPLIFICATION","DELETION"))

  expect_no_error(cna_wide <- pivot_cna_wider(cna_long))

  nrow_long <- cna_long %>%
    filter(sampleId == "P-0001276-T01-IM3" & alteration == "AMPLIFICATION") %>%
    nrow(.)

  len_wide <- length(cna_wide$`P-0001276-T01-IM3`[cna_wide$`P-0001276-T01-IM3` == 2])

  expect_equal(nrow_long, len_wide)
})

test_that("pivot CNA- check ampl/del vs with 2/-1 alteration coding", {

  cna_long_1 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c("AMPLIFICATION","AMPLIFICATION",
                   "DELETION","AMPLIFICATION","DELETION"))

  cna_wide_1 <- pivot_cna_wider(cna_long_1)

   cna_long_2 <- data.frame(
      sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                   "P-0005436-T01-IM3",
                   "P-0001276-T01-IM3","P-0003333-T01-IM3"),
      Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                      "HIST1H3B","KDR"),
      alteration = c(2, 2, -2, 2, -2))

  cna_wide_2 <- pivot_cna_wider(cna_long_2)
  expect_equal(cna_wide_1, cna_wide_2)

  cna_long_3 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c("amplification","AMPLIFICATION",
                   "deletion","AMPLIFICATION","DELETION"))

  cna_wide_3 <- pivot_cna_wider(cna_long_3)
  expect_equal(cna_wide_3, cna_wide_2)

})

# INSERT TEST HERE
# Check above but with -1/1 once we get that final coding!!


test_that("pivot CNA- check when NAs", {

  cna_long_1 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(NA, "AMPLIFICATION",
                   "DELETION","AMPLIFICATION","DELETION"))

  cna_wide_1 <- pivot_cna_wider(cna_long_1)

  cna_long_2 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(NA, 2, -2, 2, -2))

  cna_wide_2 <- pivot_cna_wider(cna_long_2)

  expect_equal(cna_wide_1, cna_wide_2)

  expect_equal(sum(is.na(cna_wide_2$`P-0001276-T01-IM3`)), 1)
  expect_equal(sum(is.na(cna_wide_1$`P-0001276-T01-IM3`)), 1)
})

test_that("pivot CNA-  unknown alteration vaues", {


  cna_long_1 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(NA, "THING",
                   "OTHER_THING","AMPLIFICATION","DELETION"))

  expect_error(pivot_cna_wider(cna_long_1), "Unknown values in alteration field:*")



})


test_that("pivot CNA-  mixed coding in alteration column and 0/neutral", {


  cna_long_1 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(NA, 2, 0, "AMPLIFICATION","DELETION"))

  wide_cna_1 <- pivot_cna_wider(cna_long_1)

  cna_long_2 <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(NA, "AMPLIFICATION", "NEUTRAL", "AMPLIFICATION","DELETION"))

  wide_cna_2 <- pivot_cna_wider(cna_long_2)

  expect_equal(wide_cna_1, wide_cna_2)

})
`

test_that("pivot CNA- non standard sample IDs with periods", {


  cna_long_1 <- data.frame(
    sampleId = c("P.0001276.T01.IM3","P.0001276.T01.IM3",
                 "P.0005436.T01.IM3",
                 "P.0001276.T01.IM3","P.0003333.T01.IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(NA, 2, 2, 1, -1))

  expect_no_error(wide_cna_1 <- pivot_cna_wider(cna_long_1))

})

# Pivot Longer -----------------------------------------------------------------
test_that("pivot CNA- no CNA events error", {

  cna_wide <- tibble::tribble(
    ~Hugo_Symbol, ~P.0070637.T01.IM7, ~P.0042589.T01.IM6, ~P.0026544.T01.IM6, ~P.0032011.T01.IM6,
    "CRKL",                 0L,                 0L,                 0L,                 0L,
    "SCG5",                 0L,                 0L,                 0L,                 0L,
    "STK11",                 0L,                 0L,                 0L,                 0L,
    "MEN1",                 0L,                 0L,                 0L,                 0L,
    "B2M",                 0L,                 0L,                 0L,                 0L,
    "TAP2",                 0L,                 0L,                 0L,                 0L,
    "PMAIP1",                 0L,                 0L,                 0L,                 0L,
    "H3-3A",                 0L,                 0L,                 0L,                 0L,
    "H3-3B",                 0L,                 0L,                 0L,                 0L,
    "CDC73",                 0L,                 0L,                 0L,                 0L,
    "PIK3CA",                 0L,                 0L,                 0L,                 0L
  )

  expect_error(long_cna_1 <- pivot_cna_longer(cna_wide), "There are*")

})

test_that("pivot CNA- test clean_sample_ids = TRUE", {

  cna_wide <- tibble::tribble(
    ~Hugo_Symbol, ~P.0070637.T01.IM7, ~P.0042589.T01.IM6, ~P.0026544.T01.IM6, ~P.0032011.T01.IM6,
    "CRKL",                 0L,                 0L,                 0L,                 -2L,
    "SCG5",                 0L,                 0L,                 0L,                 0L,
    "STK11",                 1L,                 0L,                 0L,                 0L,
    "MEN1",                 0L,                 0L,                 0L,                 0L,
    "B2M",                 0L,                 2L,                 0L,                 0L,
    "TAP2",                 0L,                 0L,                 -1L,                 0L,
    "PMAIP1",                 0L,                 0L,                 0L,                 0L,
    "H3-3A",                 0L,                 0L,                 0L,                 0L,
    "H3-3B",                 0L,                 0L,                 0L,                 0L,
    "CDC73",                 0L,                 0L,                 0L,                 0L,
    "PIK3CA",                 0L,                 0L,                 0L,                 -1L
  )

  expect_message(long_cna_1 <- pivot_cna_longer(cna_wide), "Replacing all*")
  expect_false(any(str_detect(long_cna_1$sample_id, fixed("."))))

})

test_that("pivot CNA- test clean_sample_ids = FALSE", {

  cna_wide <- tibble::tribble(
    ~Hugo_Symbol, ~P.0070637.T01.IM7, ~P.0042589.T01.IM6, ~P.0026544.T01.IM6, ~P.0032011.T01.IM6,
    "CRKL",                 0L,                 0L,                 0L,                 -2L,
    "SCG5",                 0L,                 0L,                 0L,                 0L,
    "STK11",                 1L,                 0L,                 0L,                 0L,
    "MEN1",                 0L,                 0L,                 0L,                 0L,
    "B2M",                 0L,                 2L,                 0L,                 0L,
    "TAP2",                 0L,                 0L,                 -1L,                 0L,
    "PMAIP1",                 0L,                 0L,                 0L,                 0L,
    "H3-3A",                 0L,                 0L,                 0L,                 0L,
    "H3-3B",                 0L,                 0L,                 0L,                 0L,
    "CDC73",                 0L,                 0L,                 0L,                 0L,
    "PIK3CA",                 0L,                 0L,                 0L,                 -1L
  )

  expect_no_message(long_cna_1 <- pivot_cna_longer(cna_wide, clean_sample_ids = FALSE))
  expect_true(any(str_detect(long_cna_1$sample_id, fixed("."))))

})

test_that("pivot CNA- test clean_sample_ids = FALSE", {

  cna_wide <- tibble::tribble(
    ~Hugo_Symbol, ~P.0070637.T01.IM7, ~P.0042589.T01.IM6, ~P.0026544.T01.IM6, ~P.0032011.T01.IM6,
    "CRKL",                 0L,                 0L,                 0L,                 -2L,
    "SCG5",                 0L,                 0L,                 0L,                 0L,
    "STK11",                 1L,                 0L,                 0L,                 0L,
    "MEN1",                 0L,                 0L,                 0L,                 0L,
    "B2M",                 0L,                 2L,                 0L,                 0L,
    "TAP2",                 0L,                 0L,                 -1L,                 0L,
    "PMAIP1",                 0L,                 0L,                 0L,                 0L,
    "H3-3A",                 0L,                 0L,                 0L,                 0L,
    "H3-3B",                 0L,                 0L,                 0L,                 0L,
    "CDC73",                 0L,                 0L,                 0L,                 0L,
    "PIK3CA",                 0L,                 0L,                 0L,                 -1L
  )

  expect_no_message(long_cna_1 <- pivot_cna_longer(cna_wide, clean_sample_ids = FALSE))
  expect_true(any(str_detect(long_cna_1$sample_id, fixed("."))))

})
