library(cbioportalr)
library(tidyverse)

# This script creates and saves the IMPACT Gene Info Data Frame within the package
# The script runs as is, but is commented out in because it requires calling the cbioportal API

# Get IMPACT genes --------------------------------------------------------
# load(here::here("data", "ti_341.rda"))
# load(here::here("data", "ti_410.rda"))
# load(here::here("data", "ti_468.rda"))
#
# l <- list("ti_341" = ti_341,
#      "ti_410" = ti_410,
#      "ti_468" = ti_468)
#
# # function to extract gene name
# extract_genes <- function(impact_plat, name) {
#
#   impact_plat %>%
#     transmute(gene = str_extract(.data$V5, "[^_]+"),
#            platform = name) %>%
#     filter(gene != "Tiling") %>%
#     distinct()
# }
#
# # create data frame of genes
# impact_genes <- map2_df(l, names(l),
#                      ~extract_genes(impact_plat = .x,
#                                     name = .y))
#
# impact_genes_wide <- impact_genes %>%
#   mutate(value = "included") %>%
#   pivot_wider(names_from = platform,
#   values_from = value, values_fill = "not_included")
#
#
# # Get Aliases  --------------------------------------------------------
#
# # a function to get aliases from cbioportal API
# get_alias <- function(hugo_symbol) {
#   url_path = paste0("genes/", hugo_symbol, "/aliases")
#   res <- cbioportalr::cbp_api(url_path)
#
#   res <- res$content
#   unlist(res)
# }
#
# # get all aliases for impact genes
# get_cbioportal_db("msk_impact")
#
# all_alias <- impact_genes_wide %>%
#   select(gene) %>%
#   dplyr::mutate(alias = purrr::map(gene,
#                 ~tryCatch(get_alias(.x),
#                           error = function(e) NA_character_)))
#
#
# all_alias_table <- all_alias %>%
#   unnest(cols = c(alias))
#
#
# # gene column is the accepted name, and alias is the name we will detect an replace with accepted name
# # therefore we need to make sure there are no accepted genes in the alias column
# all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()
#
# #"FP"     "RAD54L" "BIM"    "ROS"    "HGF"    "MLL3"   "MLL"    "KRAS"   "MLL2"   "TP53"   "EZH1"
#
# # function to resolve aliases
# recode_alias <- function(data = all_alias_table, accepted_gene, alias_gene) {
#   data <- data %>%
#     mutate(gene =
#            case_when(
#              gene == alias_gene ~ accepted_gene,
#              TRUE ~ gene))
#
#   add_on <- tribble(
#     ~gene, ~alias,
#     accepted_gene, alias_gene)
#
#   data <- data %>%
#     bind_rows(., add_on)
#
#
#   return(data)
# }
#
#
# all_alias_table <- all_alias_table %>%
#   recode_alias(., accepted_gene = "MLL", alias_gene = "KMT2A") %>%
#   recode_alias(., accepted_gene = "SDHA", alias_gene = "FP") %>%
#   recode_alias(., accepted_gene = "ATRX", alias_gene = "RAD54") %>%
#   recode_alias(., accepted_gene = "BCL2L11", alias_gene ="BIM") %>%
#   recode_alias(., accepted_gene = "ROS1", alias_gene = "ROS") %>%
#   recode_alias(., accepted_gene = "MLL3", alias_gene = "KMT2C") %>%
#   recode_alias(., accepted_gene = "MLL2", alias_gene = "KMT2D") %>%
#   recode_alias(., accepted_gene = "MLL4", alias_gene = "KMT2B") %>%
#
#   # I don't think these are aliases even though they came up.
#   # I think they are separate genes although I am not sure
#   filter(!(gene == "SOS1" & alias == "HGF")) %>%
#   filter(!(gene == "HRAS" & alias == "KRAS")) %>%
#   filter(!(gene == "TP53BP1" & alias == "TP53")) %>%
#   filter(!(gene == "EZH2" & alias == "EZH1")) %>%
#
#   # are MLL2 and MLL4 interchangeable?!
#   filter(!(gene == "MLL2" & alias == "MLL4")) %>%
#   filter(!(gene == "MLL4" & alias == "MLL2")) %>%
#   filter(!(gene == "ATRX" & alias == "RAD54L")) %>%
#   filter(!(gene == alias))
#
# # check we've recoded all
# all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()
#
# all_alias_table <- all_alias_table %>%
#   distinct() %>%
#   rename("hugo_symbol" = gene)
#
# impact_genes_wide$gene <- map_chr(impact_genes_wide$gene,
#                                   ~resolve_alias(.x,
#                                                  alias_table = all_alias_table))
#
# impact_genes_wide <- impact_genes_wide %>%
#   distinct()
#
# # there are some discrepancies with included/not included but all dupes are included in all panels so fixing
# impact_genes_wide <- impact_genes_wide %>%
#   group_by(gene) %>%
#   mutate(sum = n()) %>%
#   mutate(across(contains("ti_"), ~case_when(sum > 1 ~ "included", TRUE ~ .x))) %>%
#   select(-sum) %>% distinct()
#
# # create nested version of alias table
# all_alias_table_nest <- all_alias_table %>%
#   group_by(gene) %>%
#   nest()
#
# # join impact gene dataframe and their aliases
# impact_genes_wide <- impact_genes_wide %>%
#   left_join(all_alias_table_nest) %>%
#   rename("alias" = data)
#
# impact_genes <- impact_genes_wide
#
# names(impact_genes) <- c("hugo_symbol" ,
#                           "platform_341",
#                          "platform_410",
#                          "platform_468",
#                          "alias")
#
# # check if against the old one. There are some discrepancies
# #setdiff(gnomeR::impact_gene_info$hugo_symbol, impact_genes$gene)
#
# usethis::use_data(impact_genes, overwrite = TRUE)
