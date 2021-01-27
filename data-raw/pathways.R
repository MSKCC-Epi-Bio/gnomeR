library(tidyverse)

pathways <- read_csv(here::here("data-raw", "pathways.csv")) %>%
  select(-source)

pathways$genes <- map(pathways$genes, ~unlist(strsplit(.x, split=", ")))

usethis::use_data(pathways, overwrite = TRUE)
