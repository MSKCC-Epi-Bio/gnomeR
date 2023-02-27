library(genieBPC)
library(dplyr)
library(readr)

# set_synapse_credentials()

# sv_long <- pull_data_synapse(cohort = "CRC",
#                              version = "v2.0-public")$CRC_v2.0$fusions %>%
#  filter(Tumor_Sample_Barcode %in% unique(sv_long$Tumor_Sample_Barcode)[1:50]) %>%
#  gnomeR::rename_columns()

sv_long <- readr::read_csv(here::here('data-raw', 'sv_long.csv'))

usethis::use_data(sv_long, overwrite = TRUE)
