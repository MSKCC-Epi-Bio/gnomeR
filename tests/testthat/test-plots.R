
# Test Plotting Functions-----------------------------------------------------------


# maf_viz --------
maf_viz

test_that("Plot layers match expectations",{
    samples <- as.character(unique(gnomeR::mut$Tumor_Sample_Barcode))[1:200]
    all_plots <- maf_viz(maf = mut %>% filter(Tumor_Sample_Barcode %in% samples))
    expect_equal(names(all_plots), c("varclass", "vartype", "snvclass", "samplevar", "topgenes"))

  })


