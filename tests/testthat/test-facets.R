
# test facets_dat() ----

test_that("test", {
   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:1000]

   samples.seg <- clin.sample %>%
     filter(Sample.Identifier %in% samples,
            as.numeric(as.character(Tumor.Purity)) > 30) %>%
     pull(Sample.Identifier)

   expect_error(facets_heatmap(seg = gnomeR::seg,
                           samples=samples.seg[0:10]), NA)

})

# test facets_heatmap() ----

test_that("test", {
  samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:1000]

  samples.seg <- clin.sample %>%
    filter(Sample.Identifier %in% samples,
           as.numeric(as.character(Tumor.Purity)) > 30) %>%
    pull(Sample.Identifier)

  expect_error(facets_heatmap(seg = gnomeR::seg,
                              samples=samples.seg[0:10]), NA)

  test_pl <- facets_heatmap(seg = gnomeR::seg, samples=samples.seg[0:20])
  expect_equal(class(test_pl$p), "trellis")
})

