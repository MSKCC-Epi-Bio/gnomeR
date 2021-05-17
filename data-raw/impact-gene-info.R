# library(cbioportalR)
# library(tidyverse)
# library(gnomeR)
#
# # This script creates and saves the IMPACT Gene Info Data Frame within the package
#
# # Functions --------------------------------------------------------------------
#
# # function to get gene panel data from cbioportal API
# get_panel <- function(panel_id) {
#
#   url_path <- paste0("gene-panels/", panel_id)
#
#   res <- cbp_api(url_path)
#   df <- tibble::as_tibble(res$content) %>%
#     mutate(data = map(genes, ~as_tibble(.x))) %>%
#     select(data) %>%
#     unnest(cols = data)
#
#   return(df)
# }
#
# # function clean gene names
# clean_genes <- function(impact_plat, name) {
#
#   impact_plat %>%
#     mutate(gene = str_replace_all(gene, "-", "."),
#               platform = name) %>%
#     filter(gene != "Tiling") %>%
#     distinct()
# }
#
# # a function to get aliases from cbioportal API
# get_alias <- function(hugo_symbol) {
#   url_path = paste0("genes/", hugo_symbol, "/aliases")
#   res <- cbioportalR::cbp_api(url_path)
#
#   res <- res$content
#   unlist(res)
# }
#
#
# # Get IMPACT genes -----------------------------------------------------------
#
# # *JJ files ---
# jj_341 <- read.delim(here::here("data-raw", "data_gene_panel_impact341.txt"),
#                                 header=F, skip=3) %>%
#   pivot_longer(everything()) %>%
#  transmute(gene = str_remove_all(value, "gene_list: "))
#
# jj_410 <- read.delim(here::here("data-raw", "data_gene_panel_impact410.txt"),
#                                 header=F, skip=2) %>%
#   pivot_longer(everything()) %>%
#  transmute(gene = str_remove_all(value, "gene_list: "))
#
#
# jj_468 <- read.delim(here::here("data-raw", "data_gene_panel_impact468.txt"),
#                                 header=F, skip=3) %>%
#   pivot_longer(everything()) %>%
#  transmute(gene = str_remove_all(value, "gene_list: "))
#
# l <- list("jj_341" = jj_341,
#           "jj_410" = jj_410,
#           "jj_468" = jj_468)
#
#
# # create data frame of genes
# jj_impact_genes <- map2_df(l, names(l),
#                         ~clean_genes(impact_plat = .x,
#                                        name = .y))
#
#
# # *API Panel Files ---
# get_cbioportal_db("cbioportal.mskcc.org")
#
# ti_341 <- get_panel(panel_id = "IMPACT341") %>%
#   rename(gene = hugoGeneSymbol,
#          entrez_id = entrezGeneId)
#
# ti_410 <- get_panel(panel_id = "IMPACT410") %>%
#   rename(gene = hugoGeneSymbol,
#          entrez_id = entrezGeneId)
#
# ti_468 <- get_panel(panel_id = "IMPACT468") %>%
#   rename(gene = hugoGeneSymbol,
#          entrez_id = entrezGeneId)
#
# l_api <- list("ti_341" = ti_341,
#           "ti_410" = ti_410,
#           "ti_468" = ti_468)
#
# # create data frame of genes
# impact_genes <- map2_df(l_api, names(l_api),
#                         ~clean_genes(impact_plat = .x,
#                                        name = .y))
#
# impact_genes_wide <- impact_genes %>%
#   mutate(value = "included") %>%
#   pivot_wider(names_from = platform,
#               values_from = value, values_fill = "not_included")
#
#
# # Get Gene Entrez IDs  --------------------------------------------------------
#
# # API call to get a list of genes and entrez ids
# all_genes <- get_genes()
#
# all_genes <- all_genes %>%
#   filter(type == "protein-coding") %>%
#   janitor::clean_names() %>%
#   rename("hugo_symbol" = hugo_gene_symbol)
#
# # Get Aliases  --------------------------------------------------------
#
# # get all aliases for impact genes
# all_alias <- impact_genes_wide %>%
#   select(gene) %>%
#   dplyr::mutate(alias = purrr::map(gene,
#                                    ~tryCatch(get_alias(.x),
#                                              error = function(e) NA_character_)))
# # unnested long version
# all_alias_table <- all_alias %>%
#   unnest(cols = c(alias),
#          keep_empty = TRUE)
#
#
# # compare genes in JJ's github files to API panel files
# diff <- setdiff(jj_impact_genes$gene, all_alias_table$gene)
#
# diff
# # [1] "FAM123B" "FAM175A" "FAM46C"  "MLL"     "MLL2"    "MLL3"    "MRE11A"  "MYCL1"   "PAK7"    "PARK2"   "RFWD2"
# # [12] "TCEB1"   "FAM58A"  "GTF2I"   "MLL4"    "SETD8"   "WHSC1"   "WHSC1L1"
#
# # check if genes that are not in ours (`diff) are available as aliases to other genes
# # result is 1 -"GTF2I"
# setdiff(diff, all_alias_table$alias)
#
#
# # function to resolve aliases
# recode_alias <- function(data = all_alias_table, accepted_gene, alias_gene) {
#   data <- data %>%
#     mutate(gene =
#              case_when(
#                gene == alias_gene ~ accepted_gene,
#                TRUE ~ gene))
#
#   add_on <- tribble(
#     ~gene, ~alias,
#     accepted_gene, alias_gene)
#
#   data <- data %>%
#     bind_rows(., add_on) %>%
#     filter(gene != alias | is.na(alias))
#
#
#   return(data)
# }
#
# # NOTE: On 2021-05-17 we decided to use API genes as accepted over JJ's list based on genecards.org. NO recoding here anymore
#
# # Check for genes duplicated in both accepted and alias columns
# # all_alias_table %>%
# #   filter(gene == alias) %>% View()
#
# all_alias_table <- all_alias_table %>%
#   filter(gene != alias)
#
# # check genes that are in both main gene col and alias col
# all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()
# #[1] "HGF"    "KRAS"   "RAD54L" "TP53"   "EZH1"
#
# all_alias_table <- all_alias_table %>%
#
#   # I don't think these are aliases even though they came up. They usually only come up one way
#   # I think they are separate genes although I am not sure
#   filter(!(gene == "SOS1" & alias == "HGF")) %>%
#   filter(!(gene == "HRAS" & alias == "KRAS")) %>%
#   filter(!(gene == "TP53BP1" & alias == "TP53")) %>%
#   filter(!(gene == "EZH2" & alias == "EZH1")) %>%
#   # are MLL2 and MLL4 interchangeable?!
#   filter(!(gene == "MLL2" & alias == "MLL4")) %>%
#   filter(!(gene == "MLL4" & alias == "MLL2")) %>%
#   filter(!(gene == "ATRX" & alias == "RAD54L")) %>%
#   distinct()
#
# # check we've recoded all
# all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()
#
# all_alias_table <- all_alias_table %>%
#   rename("hugo_symbol" = gene)
#
#
# impact_genes_wide <- impact_genes_wide %>%
#   distinct()
#
#
# # create nested version of alias table to be joined to main dataframe
# all_alias_table_nest <- all_alias_table %>%
#   left_join(., select(impact_genes, gene, entrez_id),
#             by = c("hugo_symbol" = "gene")) %>%
#   left_join(., select(all_genes, hugo_symbol, entrez_gene_id),
#             by = c("alias" = "hugo_symbol")) %>%
#   rename(alias_entrez_id = entrez_gene_id) %>%
#
#   group_by(hugo_symbol) %>%
#   nest(data = c(alias, alias_entrez_id))
#
#
# # join impact gene dataframe and their aliases
# impact_genes_wide <- impact_genes_wide %>%
#   rename("hugo_symbol" = gene) %>%
#   left_join(., all_alias_table_nest, by = c("hugo_symbol", "entrez_id"))
#
# impact_gene_info <- impact_genes_wide
#
# names(impact_gene_info) <- c("hugo_symbol" ,
#                          "platform_341",
#                          "platform_410",
#                          "platform_468",
#                          "entrez_id",
#                          "alias")
#
# usethis::use_data(impact_gene_info , overwrite = TRUE)

# save(impact_gene_info, file = here::here("impact_gene_info.RData"))


