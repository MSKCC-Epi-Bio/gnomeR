library(tidyverse)

pathways <- read_csv(here::here("data-raw", "pathways.csv"))

pathways <- map(pathways, ~ unique(.x[!is.na(.x)]))

usethis::use_data(pathways, overwrite = TRUE)
