library(cbioportalR)
library(dplyr)
library(readr)

names_df <- read_csv(here::here("data-raw", "accepted-column-names.csv"),
                     trim_ws = TRUE)

names_df <- names_df %>%
  mutate(sc_maf_column_name = snakecase::to_snake_case(maf_column_name)) %>%
  mutate(sc_api_column_name = snakecase::to_snake_case(api_column_name)) %>%

  mutate(caps_maf_column_name = toupper(sc_maf_column_name)) %>%
  mutate(caps_api_column_name = toupper(sc_api_column_name)) %>%

  mutate(internal_column_name = case_when(
    sc_maf_column_name == "tumor_sample_barcode" ~ "sample_id",
    sc_maf_column_name == "class" ~ "variant_class",
    is.na(sc_maf_column_name) & !is.na(sc_api_column_name) ~ sc_api_column_name,
    TRUE ~ sc_maf_column_name
  ))


usethis::use_data(names_df, overwrite = TRUE)
