library(cbioportalR)
library(dplyr)

# Public
cbioportalR::set_cbioportal_db("public")

all_public_panels <- available_gene_panels()
gene_panels <- cbioportalR::get_gene_panel(all_public_panels$genePanelId)

gene_panels <- gene_panels %>%
  transmute(gene_panel = .data$genePanelId,
            entrez_id = .data$entrezGeneId,
            hugo_symbol = .data$hugoGeneSymbol)

gene_panels <- gene_panels %>%
  group_by(.data$gene_panel) %>%
  summarise(genes_in_panel = list(.data$hugo_symbol),
            entrez_ids_in_panel = list(.data$entrez_id))



# MSK Additional ------------------------------------------------------------
# cbioportalR::set_cbioportal_db("msk")
#
# all_msk_panels <- available_gene_panels()
# additional_msk_panels <- setdiff(all_msk_panels$genePanelId, gene_panels$gene_panel)
#
# add_msk_panels <- cbioportalR::get_gene_panel(additional_msk_panels)
#
# add_msk_panels <- add_msk_panels %>%
#   transmute(gene_panel = .data$genePanelId,
#             entrez_id = .data$entrezGeneId,
#             hugo_symbol = .data$hugoGeneSymbol)
#
#
#
# add_msk_panels <- add_msk_panels %>%
#   group_by(.data$gene_panel) %>%
#   summarise(genes_in_panel = list(.data$hugo_symbol),
#             entrez_ids_in_panel = list(.data$entrez_id))
#
#
# gene_panels <- bind_rows(gene_panels, add_msk_panels)

# GENIE Additional ------------------------------------------------------------

cbioportalR::set_cbioportal_db("https://genie-private.cbioportal.org/api")

all_genie_panels <- available_gene_panels()
additional_genie_panels <- setdiff(all_genie_panels$genePanelId, gene_panels$gene_panel)

add_genie_panels <- cbioportalR::get_gene_panel(additional_genie_panels)

# get table of gene aliases (internal data)
alias_table <- gnomeR::impact_alias_table %>%
  dplyr::select("hugo_symbol", "alias")

# recode aliases
add_genie_panels$hugoGeneSymbol <- purrr::map_chr(add_genie_panels$hugoGeneSymbol,
                                         ~resolve_alias(gene_to_check = .x,
                                                        alias_table = alias_table))

add_genie_panels <- add_genie_panels %>%
  transmute(gene_panel = .data$genePanelId,
            entrez_id = .data$entrezGeneId,
            hugo_symbol = .data$hugoGeneSymbol)



add_genie_panels <- add_genie_panels %>%
  group_by(.data$gene_panel) %>%
  summarise(genes_in_panel = list(.data$hugo_symbol),
            entrez_ids_in_panel = list(.data$entrez_id))


gene_panels <- bind_rows(gene_panels, add_genie_panels)

usethis::use_data(gene_panels, overwrite = TRUE)

