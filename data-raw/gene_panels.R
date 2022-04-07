library(cbioportalR)
library(dplyr)

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


usethis::use_data(gene_panels, overwrite = TRUE)
