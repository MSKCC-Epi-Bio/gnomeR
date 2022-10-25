library(cbioportalR)
library(dplyr)

set_cbioportal_db("public")
gen <- get_genetics_by_study("prad_broad_2013")

cbp_mut <- gen$mutation
cbp_cna <- gen$cna
cbp_sv <- gen$structural_variant

usethis::use_data(cbp_mut, overwrite = TRUE)
usethis::use_data(cbp_cna, overwrite = TRUE)
usethis::use_data(cbp_sv, overwrite = TRUE)
