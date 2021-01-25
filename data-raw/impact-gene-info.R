library(cbioportalr)
library(tidyverse)

# This script creates and saves the IMPACT Gene Info Data Frame within the package
# The script runs as is, but is commented out in because it requires calling the cbioportal API

# Get IMPACT genes --------------------------------------------------------
ti_341 <- read.delim(here::here("data-raw", "data_gene_panel_impact341.txt"),
                                header=F, skip=3) %>%
  pivot_longer(everything()) %>%
 transmute(gene = str_remove_all(value, "gene_list: "))

ti_410 <- read.delim(here::here("data-raw", "data_gene_panel_impact410.txt"),
                                header=F, skip=2) %>%
  pivot_longer(everything()) %>%
 transmute(gene = str_remove_all(value, "gene_list: "))


ti_468 <- read.delim(here::here("data-raw", "data_gene_panel_impact468.txt"),
                                header=F, skip=3) %>%
  pivot_longer(everything()) %>%
 transmute(gene = str_remove_all(value, "gene_list: "))



l <- list("ti_341" = ti_341,
          "ti_410" = ti_410,
          "ti_468" = ti_468)

# function to extract gene name
clean_genes <- function(impact_plat, name) {

  impact_plat %>%
    transmute(gene = str_replace_all(gene, "-", "."),
              platform = name) %>%
    filter(gene != "Tiling") %>%
    distinct()
}

# create data frame of genes
impact_genes <- map2_df(l, names(l),
                        ~clean_genes(impact_plat = .x,
                                       name = .y))

impact_genes_wide <- impact_genes %>%
  mutate(value = "included") %>%
  pivot_wider(names_from = platform,
              values_from = value, values_fill = "not_included")


# manually add on and recode MLL* to  KMT2* family because currently in cbioportal aliases only come up for KMT2* family and not for MLL*
# This seems counter-intuitive to code it wrong but this is so the alias function can recognize any aliases for this gene- no aliases come up for MLL* family (yet- they will likely change this)
impact_genes_wide  <- impact_genes_wide %>%
  mutate(gene = case_when(
    gene ==  "MLL" ~ "KMT2A",
    gene == "MLL2" ~ "KMT2B",
    gene == "MLL3" ~  "KMT2C",
    gene == "MLL4" ~ "KMT2D",
    TRUE ~ gene
    ))

# Get Aliases  --------------------------------------------------------

# a function to get aliases from cbioportal API
get_alias <- function(hugo_symbol) {
  url_path = paste0("genes/", hugo_symbol, "/aliases")
  res <- cbioportalr::cbp_api(url_path)

  res <- res$content
  unlist(res)
}

# get all aliases for impact genes
get_cbioportal_db("msk_impact")

all_alias <- impact_genes_wide %>%
  select(gene) %>%
  dplyr::mutate(alias = purrr::map(gene,
                                   ~tryCatch(get_alias(.x),
                                             error = function(e) NA_character_)))


all_alias_table <- all_alias %>%
  unnest(cols = c(alias))

# add on MLL
all_alias_table <- all_alias_table  %>%
  bind_rows(
    tribble(
      ~gene, ~alias,
      "MLL", "KMT2A",
      "MLL2", "KMT2B",
      "MLL3", "KMT2C",
      "MLL4", "KMT2D"
      ))


# gene column is the accepted name, and alias is the name we will detect an replace with accepted name
# therefore we need to make sure there are no accepted genes in the alias column
all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()

# "FAM123B" "HGF"     "KRAS"    "MLL"     "MLL2"    "MLL3"    "RAD54L"  "TP53"    "EZH1"    "MLL4"

# function to resolve aliases
recode_alias <- function(data = all_alias_table, accepted_gene, alias_gene) {
  data <- data %>%
    mutate(gene =
             case_when(
               gene == alias_gene ~ accepted_gene,
               TRUE ~ gene))

  add_on <- tribble(
    ~gene, ~alias,
    accepted_gene, alias_gene)

  data <- data %>%
    bind_rows(., add_on)


  return(data)
}


all_alias_table <- all_alias_table %>%
  recode_alias(., accepted_gene = "FAM123B", alias_gene = "AMER1") %>%
  recode_alias(., accepted_gene = "MLL", alias_gene = "KMT2A") %>%
  recode_alias(., accepted_gene = "MLL2", alias_gene = "KMT2B") %>%
  recode_alias(., accepted_gene = "MLL3", alias_gene = "KMT2C") %>%
  recode_alias(., accepted_gene = "MLL4", alias_gene = "KMT2D") %>%

  # I don't think these are aliases even though they came up.
  # I think they are separate genes although I am not sure
  filter(!(gene == "SOS1" & alias == "HGF")) %>%
  filter(!(gene == "HRAS" & alias == "KRAS")) %>%
  filter(!(gene == "TP53BP1" & alias == "TP53")) %>%
  filter(!(gene == "EZH2" & alias == "EZH1")) %>%
  # are MLL2 and MLL4 interchangeable?!
  filter(!(gene == "MLL2" & alias == "MLL4")) %>%
  filter(!(gene == "MLL4" & alias == "MLL2")) %>%
  filter(!(gene == "ATRX" & alias == "RAD54L")) %>%
  filter(!(gene == alias)) %>%

  distinct()



# check we've recoded all
all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()

all_alias_table <- all_alias_table %>%
  distinct() %>%
  rename("hugo_symbol" = gene)

impact_genes_wide$gene <- map_chr(impact_genes_wide$gene,
                                  ~resolve_alias(.x,
                                                 alias_table = all_alias_table))

impact_genes_wide <- impact_genes_wide %>%
  distinct()

# there are some discrepancies with included/not included but all dupes are included in all panels so fixing
impact_genes_wide <- impact_genes_wide %>%
  group_by(gene) %>%
  mutate(sum = n()) %>%
  mutate(across(contains("ti_"), ~case_when(sum > 1 ~ "included", TRUE ~ .x))) %>%
  select(-sum) %>% distinct()

# create nested version of alias table
all_alias_table_nest <- all_alias_table %>%
  group_by(hugo_symbol) %>%
  nest()

# join impact gene dataframe and their aliases
impact_genes_wide <- impact_genes_wide %>%
  rename(hugo_symbol = gene) %>%
  left_join(all_alias_table_nest) %>%
  rename("alias" = data)

impact_gene_info <- impact_genes_wide

names(impact_gene_info) <- c("hugo_symbol" ,
                         "platform_341",
                         "platform_410",
                         "platform_468",
                         "alias")

usethis::use_data(impact_gene_info , overwrite = TRUE)
