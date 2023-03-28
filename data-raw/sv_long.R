#library(genieBPC)
library(dplyr)
library(readr)

# set_synapse_credentials()

# sv_long <- pull_data_synapse(cohort = "NSCLC",
#                              version = "v2.0-public")$NSCLC_v2.0$fusions
#
# samples <- sample(x = sv_long$Tumor_Sample_Barcode, size = 30, replace = F)
#
# sv_long <- sv_long %>%
#  filter(Tumor_Sample_Barcode %in% samples) %>%
#  gnomeR::rename_columns()
#
# readr::write_csv(sv_long, here::here('data-raw', 'sv_long.csv'))


sv_long <- readr::read_csv(here::here('data-raw', 'sv_long.csv'))

usethis::use_data(sv_long, overwrite = TRUE)

