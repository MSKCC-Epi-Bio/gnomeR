# Package index

### Processing Data

- [`create_gene_binary()`](https://mskcc-epi-bio.github.io/gnomeR/reference/create_gene_binary.md)
  : Enables creation of a binary matrix from a mutation, fusion or CNA
  file with a predefined list of samples (rows are samples and columns
  are genes)

- [`summarize_by_gene()`](https://mskcc-epi-bio.github.io/gnomeR/reference/summarize_by_gene.md)
  : Simplify binary matrix to one column per gene that counts any
  alteration type as 1

- [`summarize_by_patient()`](https://mskcc-epi-bio.github.io/gnomeR/reference/summarize_by_patient.md)
  : Simplify binary matrix to one column per patient that counts any
  alteration type across all samples as 1

- [`pivot_cna_wider()`](https://mskcc-epi-bio.github.io/gnomeR/reference/pivot_cna_wider.md)
  : Pivot CNA from maf (long) version to wide version

- [`pivot_cna_longer()`](https://mskcc-epi-bio.github.io/gnomeR/reference/pivot_cna_longer.md)
  : Reformat Wide CNA Data to Long

- [`add_pathways()`](https://mskcc-epi-bio.github.io/gnomeR/reference/add_pathways.md)
  : Pathway Alterations

- [`recode_alias()`](https://mskcc-epi-bio.github.io/gnomeR/reference/recode_alias.md)
  : Recode Hugo Symbol Column

- [`reformat_fusion()`](https://mskcc-epi-bio.github.io/gnomeR/reference/reformat_fusion.md)
  :

  Enables users to reformat fusions files so that each fusion is listed
  as one row with two hugo-symbol sites instead of two rows, one for
  each site. This is the required format for the `create_gene_binary`
  function.

- [`subset_by_panel()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_panel.md)
  : Subset a Binary Matrix By Genes Available on Specified Panel

### Analyzing Data

- [`tbl_genomic()`](https://mskcc-epi-bio.github.io/gnomeR/reference/tbl_genomic.md)
  : tbl_genomic
- [`mutation_viz()`](https://mskcc-epi-bio.github.io/gnomeR/reference/mutation_viz.md)
  : Creates a set of plot summarising a mutation file.
- [`ggvarclass()`](https://mskcc-epi-bio.github.io/gnomeR/reference/ggvarclass.md)
  : Barplot of Variant Classification Counts
- [`ggvartype()`](https://mskcc-epi-bio.github.io/gnomeR/reference/ggvartype.md)
  : Barplot of Variant Type Counts
- [`ggsamplevar()`](https://mskcc-epi-bio.github.io/gnomeR/reference/ggsamplevar.md)
  : \#' Utility Function to Extract SNV \#' \#' @param x string \#'
  @param n number of characters from right \#' \#' @return string \#'
  @noRd \#' @examples \#' substrRight("Hello", 2)
- [`ggtopgenes()`](https://mskcc-epi-bio.github.io/gnomeR/reference/ggtopgenes.md)
  : Barplot of Most Frequently Altered Genes
- [`gggenecor()`](https://mskcc-epi-bio.github.io/gnomeR/reference/gggenecor.md)
  : Correlation Heatmap of the Top Altered Genes
- [`ggcomut()`](https://mskcc-epi-bio.github.io/gnomeR/reference/ggcomut.md)
  : Comutation Heatmap of the Top Altered Genes
- [`subset_by_frequency()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_frequency.md)
  : Subset a Binary Matrix By Alteration Frequency Threshold

### Helper Functions

- [`which_impact_panel()`](https://mskcc-epi-bio.github.io/gnomeR/reference/which_impact_panel.md)
  : provide a list of impact panels a provided gene is found within
- [`recode_cna()`](https://mskcc-epi-bio.github.io/gnomeR/reference/recode_cna.md)
  : Function to recode numeric CNA alteration values to factor values
- [`rename_columns()`](https://mskcc-epi-bio.github.io/gnomeR/reference/rename_columns.md)
  : Rename columns from API results to work with gnomeR functions
- [`resolve_alias()`](https://mskcc-epi-bio.github.io/gnomeR/reference/resolve_alias.md)
  : Resolve Hugo Symbol Names with Aliases
- [`extract_patient_id()`](https://mskcc-epi-bio.github.io/gnomeR/reference/extract_patient_id.md)
  : Extract IMPACT Patient ID From Sample ID

### Color Palette

- [`gnomer_colors`](https://mskcc-epi-bio.github.io/gnomeR/reference/gnomer_colors.md)
  : List of suggested color palettes for when you need a large palette
- [`gnomer_palettes`](https://mskcc-epi-bio.github.io/gnomeR/reference/gnomer_palettes.md)
  : Complete list of gnomeR color palettes
- [`gnomer_palette()`](https://mskcc-epi-bio.github.io/gnomeR/reference/gnomer_palette.md)
  : Access the colors in a gnomeR color palette
- [`set_gnomer_palette()`](https://mskcc-epi-bio.github.io/gnomeR/reference/set_gnomer_palette.md)
  : Set gnomeR color palette
- [`reset_gnomer_palette()`](https://mskcc-epi-bio.github.io/gnomeR/reference/reset_gnomer_palette.md)
  : Reset gnomeR color palette

### Example Data Sets

- [`mutations`](https://mskcc-epi-bio.github.io/gnomeR/reference/mutations.md)
  : An example IMPACT cBioPortal mutation data set in API format
- [`cna`](https://mskcc-epi-bio.github.io/gnomeR/reference/cna.md) : An
  example IMPACT cBioPortal mutation data set in API format
- [`sv`](https://mskcc-epi-bio.github.io/gnomeR/reference/sv.md) : An
  example IMPACT cBioPortal mutation data set in API format
- [`seg`](https://mskcc-epi-bio.github.io/gnomeR/reference/seg.md) : A
  segmentation file from the cbioPortal datasets
- [`cna_wide`](https://mskcc-epi-bio.github.io/gnomeR/reference/cna_wide.md)
  : An example IMPACT cBioPortal CNA in wide format
- [`consequence_map`](https://mskcc-epi-bio.github.io/gnomeR/reference/consequence_map.md)
  : Consequence Map
- [`gene_panels`](https://mskcc-epi-bio.github.io/gnomeR/reference/gene_panels.md)
  : Public Gene Panels on cBioPortal
- [`names_df`](https://mskcc-epi-bio.github.io/gnomeR/reference/names_df.md)
  : Data Frame of Column Names
- [`pathways`](https://mskcc-epi-bio.github.io/gnomeR/reference/pathways.md)
  : IMPACT Oncogenic Signaling Pathways
- [`impact_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/impact_alias_table.md)
  : IMPACT Alias Tables
- [`genie_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/genie_alias_table.md)
  : GENIE Alias Table
- [`sv_long`](https://mskcc-epi-bio.github.io/gnomeR/reference/sv_long.md)
  : An example of long-format fusion/sv files
- [`clin_collab_df`](https://mskcc-epi-bio.github.io/gnomeR/reference/clin_collab_df.md)
  : An example data set for an IMPACT analysis coming from a clinical
  collaborator
- [`genie_mut`](https://mskcc-epi-bio.github.io/gnomeR/reference/genie_mut.md)
  : An example GENIE BPC mutations data set
- [`genie_cna`](https://mskcc-epi-bio.github.io/gnomeR/reference/genie_cna.md)
  : An example GENIE BPC CNA data set
- [`genie_fusion`](https://mskcc-epi-bio.github.io/gnomeR/reference/genie_fusion.md)
  : An example GENIE BPC fusions data set
