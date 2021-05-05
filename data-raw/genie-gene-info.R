# library(cbioportalR)
# library(tidyverse)
# library(gnomeR)
# library(GenieBPC)
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
#            platform = name) %>%
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
#                      header=F, skip=3) %>%
#   pivot_longer(everything()) %>%
#   transmute(gene = str_remove_all(value, "gene_list: "))
#
# jj_410 <- read.delim(here::here("data-raw", "data_gene_panel_impact410.txt"),
#                      header=F, skip=2) %>%
#   pivot_longer(everything()) %>%
#   transmute(gene = str_remove_all(value, "gene_list: "))
#
#
# jj_468 <- read.delim(here::here("data-raw", "data_gene_panel_impact468.txt"),
#                      header=F, skip=3) %>%
#   pivot_longer(everything()) %>%
#   transmute(gene = str_remove_all(value, "gene_list: "))
#
# l <- list("jj_341" = jj_341,
#           "jj_410" = jj_410,
#           "jj_468" = jj_468)
#
#
# # create data frame of genes
# jj_impact_genes <- map2_df(l, names(l),
#                            ~clean_genes(impact_plat = .x,
#                                         name = .y))
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
# l_jj <- list("ti_341" = ti_341,
#              "ti_410" = ti_410,
#              "ti_468" = ti_468)
#
# # create data frame of genes
# impact_genes <- map2_df(l_jj, names(l_jj),
#                         ~clean_genes(impact_plat = .x,
#                                      name = .y))
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
# #[1] "FAM123B" "FAM175A" "FAM46C"  "MLL"     "MLL2"    "MLL3"    "MRE11A"  "MYCL1"
# # [9] "PAK7"    "PARK2"   "RFWD2"
#
# #  [1] "FAM123B" "FAM175A" "FAM46C"  "MLL"     "MLL2"    "MLL3"    "MRE11A"  "MYCL1"
# #  [9] "PAK7"    "PARK2"   "RFWD2"
# #  "TCEB1"   "FAM58A"  "GTF2I"   "MLL4"    "SETD8"
# # [17] "WHSC1"   "WHSC1L1"
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
# # chosing the ones in JJ's list as the "main" and the API as "alias" under the
# # assumption that JJ's list is more up to date because it has MLL family
# all_alias_table <- all_alias_table %>%
#   recode_alias(., accepted_gene = "FAM123B", alias_gene = "AMER1") %>%
#   recode_alias(., accepted_gene = "FAM175A", alias_gene = "ABRAXAS1") %>%
#   recode_alias(., accepted_gene = "FAM46C", alias_gene = "TENT5C") %>%
#
#   recode_alias(., accepted_gene = "MLL", alias_gene = "KMT2A") %>%
#   recode_alias(., accepted_gene = "MLL2", alias_gene = "KMT2D") %>%
#   recode_alias(., accepted_gene = "MLL3", alias_gene = "KMT2C") %>%
#   recode_alias(., accepted_gene = "MLL4", alias_gene = "KMT2B") %>%
#
#   recode_alias(., accepted_gene = "MRE11A", alias_gene = "MRE11") %>%
#   recode_alias(., accepted_gene = "MYCL1", alias_gene = "MYCL") %>%
#   recode_alias(., accepted_gene = "PAK7", alias_gene = "PAK5") %>%
#   recode_alias(., accepted_gene = "PARK2", alias_gene = "PRKN") %>%
#   recode_alias(., accepted_gene = "RFWD2", alias_gene = "COP1") %>%
#   recode_alias(., accepted_gene = "TCEB1", alias_gene = "ELOC") %>%
#   recode_alias(., accepted_gene = "FAM58A", alias_gene = "CCNQ") %>%
#   recode_alias(., accepted_gene = "SETD8", alias_gene = "KMT5A") %>%
#   recode_alias(., accepted_gene = "WHSC1", alias_gene = "NSD2") %>%
#   recode_alias(., accepted_gene = "WHSC1L1", alias_gene = "NSD3")
#
# # check genes that are in both main gene col and alias col
# all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()
# #[1] "HGF"    "MLL4"   "KRAS"   "RAD54L" "TP53"   "EZH1"   "MLL2"
#
# all_alias_table <- all_alias_table %>%
#
#   # I don't think these are aliases even though they came up.
#   # I think they are separate genes although I am not sure
#   filter(!(gene == "SOS1" & alias == "HGF")) %>%
#   filter(!(gene == "HRAS" & alias == "KRAS")) %>%
#   filter(!(gene == "TP53BP1" & alias == "TP53")) %>%
#   filter(!(gene == "EZH2" & alias == "EZH1")) %>%
#   # are MLL2 and MLL4 interchangeable?!
#   filter(!(gene == "MLL2" & alias == "MLL4")) %>%
#   filter(!(gene == "MLL4" & alias == "MLL2")) %>%
#   filter(!(gene == "ATRX" & alias == "RAD54L")) %>%
#   filter(!(gene == alias) | is.na(alias)) %>%
#   distinct()
#
# # check we've recoded all
# all_alias_table$gene[map_lgl(all_alias_table$gene, ~.x %in% all_alias_table$alias)] %>% unique()
#
# all_alias_table <- all_alias_table %>%
#   rename("hugo_symbol" = gene)
#
# impact_genes_wide$gene <- map_chr(impact_genes_wide$gene,
#                                   ~resolve_alias(.x,
#                                                  alias_table = all_alias_table))
#
# impact_genes_wide <- impact_genes_wide %>%
#   distinct()
#
# # there are some discrepancies with included/not included in panels but all dupes are included in all panels so fixing
# impact_genes_wide <- impact_genes_wide %>%
#   group_by(gene) %>%
#   mutate(sum = n()) %>%
#   mutate(across(contains("ti_"), ~case_when(sum > 1 ~ "included", TRUE ~ .x))) %>%
#   select(-sum) %>% distinct() %>%
#
#   # need to freshly rejoin entrez ID because recoding above
#   select(-entrez_id)
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
# # if entrez id is NA, grab the first alias entrez ID (if available) so we have an ID for everything
# all_alias_table_nest  <- all_alias_table_nest %>%
#   mutate(entrez_id = case_when(
#     is.na(entrez_id) ~ map(data, ~na.omit(.x$alias_entrez_id)[1]),
#     TRUE ~ list(entrez_id)))
#
#
#
# # join impact gene dataframe and their aliases
# impact_genes_wide <- impact_genes_wide %>%
#   rename("hugo_symbol" = gene) %>%
#   left_join(., all_alias_table_nest)
#
# impact_gene_info <- impact_genes_wide
#
# names(impact_gene_info) <- c("hugo_symbol" ,
#                              "platform_341",
#                              "platform_410",
#                              "platform_468",
#                              "entrez_id",
#                              "alias")
#
#
# #########################################
#
# # adding gene panels from GENIE centers #
# # get all existing genes in the current database #
# all_gene_names <- unique(
#   c(impact_gene_info$hugo_symbol,
#     unique(as.character(unlist(impact_gene_info %>% select(alias) %>% unnest(cols = c(alias)) %>% select(alias))))
#   )
# )
# all_gene_names <- all_gene_names[!is.na(all_gene_names)]
#
# # get all genes that exists in the genie centers #
# keep_panels <- panel_names$Panel[grep("MSK",panel_names$Panel,invert = T)]
# genes <- unique(unlist(lapply(keep_panels,function(x){
#   get(x)
# })))
#
# # transform all gene names that are aliases into the accepted hugo_symbol #
# genes2 <- unique(unlist(lapply(genes, function(g){
#   if(g == "INSRR")
#     g <- "INSR"
#   if(g == "GNB2L1")
#     g <- "RACK1"
#   if(g == "PVRL4")
#     g <- "NECTIN4"
#   if(g == "C1orf86")
#     g <- "FAAP20"
#   if(g == "C17orf70")
#     g <- "FAAP100"
#   if(g == "C19orf40")
#     g <- "FAAP24"
#   if(g == "BRE")
#     g <- "BABAM2"
#   if(g == "WISP3")
#     g <- "CCN6"
#   if(g == "GPR124")
#     g <- "ADGRA2"
#   if(g == "C11orf30")
#     g <- "EMSY"
#
#   if(g %in% impact_gene_info$hugo_symbol)
#     return(g)
#   ind <- match(g, all_gene_names)
#   if(!is.na(ind)){
#     # print(ind)
#     temp_gene_alias <- as.character(unlist(impact_gene_info %>%
#                                              unnest(cols = 'alias') %>%
#                                              filter(alias == g) %>%
#                                              select(hugo_symbol)))
#
#     return(temp_gene_alias)
#   }
#   return(g)
# })))
#
# # lapply(genes2, function(gg){
# #   print(gg)
# #   get_gene_id(gg)
# # })
#
# # transform panels accordingly #
# panels <- lapply(keep_panels,function(x){
#   unique(unlist(lapply(get(x), function(g){
#     if(g == "INSRR")
#       g <- "INSR"
#     if(g == "GNB2L1")
#       g <- "RACK1"
#     if(g == "PVRL4")
#       g <- "NECTIN4"
#     if(g == "C1orf86")
#       g <- "FAAP20"
#     if(g == "C17orf70")
#       g <- "FAAP100"
#     if(g == "C19orf40")
#       g <- "FAAP24"
#     if(g == "BRE")
#       g <- "BABAM2"
#     if(g == "WISP3")
#       g <- "CCN6"
#     if(g == "GPR124")
#       g <- "ADGRA2"
#     if(g == "C11orf30")
#       g <- "EMSY"
#     if(g %in% impact_gene_info$hugo_symbol)
#       return(g)
#     ind <- match(g, all_gene_names)
#     if(!is.na(ind)){
#       # print(ind)
#       temp_gene_alias <- as.character(unlist(impact_gene_info %>%
#                                                unnest(cols = 'alias') %>%
#                                                filter(alias == g) %>%
#                                                select(hugo_symbol)))
#
#       return(temp_gene_alias)
#     }
#     return(g)
#   }
#
#   )))
# })
# names(panels) <- keep_panels
#
#
# # add all these panels to impact_gene_info #
# test_1 <- lapply(genes2, function(x){
#   print(x)
#   if(x == "INSRR")
#     x <- "INSR"
#
#   # try to find in accepted hugo symbols #
#   ind <- match(x, impact_gene_info$hugo_symbol)
#   if(!is.na(ind)){
#     temp <- impact_gene_info %>% filter(hugo_symbol == x) %>%
#       mutate(
#         DFCI_1 = ifelse(x %in% panels$DFCI_1,"included", "not_included"),
#         DFCI_2 = ifelse(x %in% panels$DFCI_2,"included", "not_included"),
#         DFCI_3 = ifelse(x %in% panels$DFCI_3,"included", "not_included"),
#         DFCI_3.1 = ifelse(x %in% panels$DFCI_3.1,"included", "not_included"),
#         UHN_48 = ifelse(x %in% panels$UHN_48,"included", "not_included"),
#         UHN_50 = ifelse(x %in% panels$UHN_50,"included", "not_included"),
#         VICC_1_SOLIDTUMOR = ifelse(x %in% panels$VICC_1_SOLIDTUMOR,"included", "not_included"),
#         VICC_1_T7 = ifelse(x %in% panels$VICC_1_T7,"included", "not_included"),
#         VICC_1_T5A = ifelse(x %in% panels$VICC_1_T5A,"included", "not_included")
#       ) %>%
#       select(hugo_symbol, platform_341, platform_410, platform_468,
#              DFCI_1, DFCI_2, DFCI_3, DFCI_3.1,
#              UHN_48, UHN_50, VICC_1_SOLIDTUMOR,
#              VICC_1_T7, VICC_1_T5A,
#              entrez_id, alias)
#     return(temp)
#   }
#
#   # try to find in aliases --> shouldn't exists but ...#
#   ind_alias <- match(x, all_gene_names)
#   if(!is.na(ind_alias)){
#     print(paste0("An alias was found:", x))
#     temp_gene_alias <- as.character(unlist(impact_gene_info %>%
#                                              unnest(cols = 'alias') %>%
#                                              filter(alias == x) %>%
#                                              select(hugo_symbol)))
#
#     temp <- impact_gene_info %>%
#       filter(hugo_symbol == temp_gene_alias) %>%
#       mutate(
#         DFCI_1 = ifelse(x %in% panels$DFCI_1,"included", "not_included"),
#         DFCI_2 = ifelse(x %in% panels$DFCI_2,"included", "not_included"),
#         DFCI_3 = ifelse(x %in% panels$DFCI_3,"included", "not_included"),
#         DFCI_3.1 = ifelse(x %in% panels$DFCI_3.1,"included", "not_included"),
#         UHN_48 = ifelse(x %in% panels$UHN_48,"included", "not_included"),
#         UHN_50 = ifelse(x %in% panels$UHN_50,"included", "not_included"),
#         VICC_1_SOLIDTUMOR = ifelse(x %in% panels$VICC_1_SOLIDTUMOR,"included", "not_included"),
#         VICC_1_T7 = ifelse(x %in% panels$VICC_1_T7,"included", "not_included"),
#         VICC_1_T5A = ifelse(x %in% panels$VICC_1_T5A,"included", "not_included")
#       ) %>%
#       select(hugo_symbol, platform_341, platform_410, platform_468,
#              DFCI_1, DFCI_2, DFCI_3, DFCI_3.1,
#              UHN_48, UHN_50, VICC_1_SOLIDTUMOR,
#              VICC_1_T7, VICC_1_T5A,
#              entrez_id, alias)
#     return(temp)
#   }
#
#   # if does not exists add new row #
#   temp_info <- get_gene_id(gsub("\\.","-",x))
#   temp <- c(x, "not_included","not_included","not_included",
#             temp_info$entrezGeneId,
#             list(get_alias(gsub("\\.","-",x))),
#             # list(get_alias(gsub("\\.","-",x))),
#             list(unlist(lapply(gsub(" ","",get_alias(gsub("\\.","-",x))),function(w){
#               alias_entrezID <- try(get_gene_id(w),silent = T)
#               if(alias_entrezID[1] %in% c(
#                 "Error : GitHub API request failed: 404\n",
#                 "Error : GitHub API request failed: 500\n"))
#                 return(NA)
#               return(alias_entrezID$entrezGeneId)
#             }))),
#             DFCI_1 = ifelse(x %in% panels$DFCI_1,"included", "not_included"),
#             DFCI_2 = ifelse(x %in% panels$DFCI_2,"included", "not_included"),
#             DFCI_3 = ifelse(x %in% panels$DFCI_3,"included", "not_included"),
#             DFCI_3.1 = ifelse(x %in% panels$DFCI_3.1,"included", "not_included"),
#             UHN_48 = ifelse(x %in% panels$UHN_48,"included", "not_included"),
#             UHN_50 = ifelse(x %in% panels$UHN_50,"included", "not_included"),
#             VICC_1_SOLIDTUMOR = ifelse(x %in% panels$VICC_1_SOLIDTUMOR,"included", "not_included"),
#             VICC_1_T7 = ifelse(x %in% panels$VICC_1_T7,"included", "not_included"),
#             VICC_1_T5A = ifelse(x %in% panels$VICC_1_T5A,"included", "not_included"))
#
#   temp <- lapply(temp,function(xx){
#     if(is.null(xx))
#       return(NA)
#     xx
#   })
#   names(temp) <- c("hugo_symbol", "platform_341", "platform_410",
#                    "platform_468", "entrez_id", "alias",
#                    "alias_entrez_id", "DFCI_1", "DFCI_2", "DFCI_3",
#                    "DFCI_3.1", "UHN_48", "UHN_50", "VICC_1_SOLIDTUMOR",
#                    "VICC_1_T7", "VICC_1_T5A")
#
#   temp <- as_tibble(temp) %>%
#     group_by(hugo_symbol) %>%
#     nest(alias = c(alias, alias_entrez_id))
#   temp <- temp %>%
#     select(hugo_symbol, platform_341, platform_410, platform_468,
#            DFCI_1, DFCI_2, DFCI_3, DFCI_3.1,
#            UHN_48, UHN_50, VICC_1_SOLIDTUMOR,
#            VICC_1_T7, VICC_1_T5A,
#            entrez_id, alias)
#   return(temp)
# })
#
# test_bind <- data.table::rbindlist(test_1)
# final <- full_join(
#   impact_gene_info,
#   test_bind,
#   by = colnames(impact_gene_info)
# ) %>%
#   mutate(DFCI_1 = ifelse(is.na(DFCI_1), "not_included",DFCI_1),
#          DFCI_2 = ifelse(is.na(DFCI_2), "not_included",DFCI_2),
#          DFCI_3 = ifelse(is.na(DFCI_3), "not_included",DFCI_3),
#          DFCI_3.1 = ifelse(is.na(DFCI_3.1), "not_included",DFCI_3.1),
#          UHN_48 = ifelse(is.na(UHN_48), "not_included",UHN_48),
#          UHN_50 = ifelse(is.na(UHN_50), "not_included",UHN_50),
#          VICC_1_SOLIDTUMOR = ifelse(is.na(VICC_1_SOLIDTUMOR), "not_included",VICC_1_SOLIDTUMOR),
#          VICC_1_T7 = ifelse(is.na(VICC_1_T7), "not_included",VICC_1_T7),
#          VICC_1_T5A = ifelse(is.na(VICC_1_T5A), "not_included",VICC_1_T5A)
#          ) %>%
#   select(hugo_symbol, platform_341, platform_410, platform_468,
#          DFCI_1,DFCI_2,DFCI_3,DFCI_3.1,UHN_48,UHN_50,
#          VICC_1_SOLIDTUMOR,VICC_1_T7,VICC_1_T5A,
#          entrez_id,alias) %>%
#   rename(MSK_341 = platform_341, MSK_410 = platform_410, MSK_468 = platform_468)
# genie_gene_info <- final
#
# usethis::use_data(genie_gene_info , overwrite = TRUE)
# usethis::use_data(panel_names , overwrite = TRUE)
#
# save(genie_gene_info, file = here::here("genie_gene_info.RData"))
# save(panel_names, file = here::here("panel_names.RData"))
#
# ########################################
#
#
#
#
