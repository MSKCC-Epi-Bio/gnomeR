
# Test Plotting Functions-----------------------------------------------------------


# maf_viz --------

test_that("Plot layers match expectations",{
    samples <- as.character(unique(gnomeR::mut$Tumor_Sample_Barcode))[1:200]

    all_plots <- mutation_viz(mutation = mut %>% filter(Tumor_Sample_Barcode %in% samples))

    expect_equal(names(all_plots), c("varclass",
                                     "vartype",
                                     "snvclass",
                                     "samplevar",
                                     "topgenes",
                                     "genecor"))
    }
)



# add.perc --------
test_that("test add_perc utility", {

  expect_error(ggvarclass(mut) + add.perc(), NA)
}

)
