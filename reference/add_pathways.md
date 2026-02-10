# Pathway Alterations

Please note that only `sample_id column`, and columns with .Amp, .Del,
.fus or no suffix are accepted. Any gene column with no suffix will be
assumed to be a mutation.

## Usage

``` r
add_pathways(
  gene_binary,
  pathways = c(names(gnomeR::pathways)),
  custom_pathways = NULL,
  other_vars = NULL,
  count_pathways_by = deprecated()
)
```

## Source

Sanchez-Vega, F., Mina, M., Armenia, J., Chatila, W. K., Luna, A., La,
K. C., Dimitriadoy, S., Liu, D. L., Kantheti, H. S., Saghafinia, S.,
Chakravarty, D., Daian, F., Gao, Q., Bailey, M. H., Liang, W. W., Foltz,
S. M., Shmulevich, I., Ding, L., Heins, Z., Ochoa, A., … Schultz, N.
(2018). Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell,
173(2), 321–337.e10. <https://doi.org/10.1016/j.cell.2018.03.035>

## Arguments

- gene_binary:

  a binary matrix from
  [`create_gene_binary()`](https://mskcc-epi-bio.github.io/gnomeR/reference/create_gene_binary.md)

- pathways:

  a vector of pre-coded pathways to annotate. The options are
  `names(gnomeR::pathways)` ("RTK/RAS", "Nrf2", "PI3K", "TGFB", "p53",
  "Wnt", "Myc", "Cell cycle", "Hippo", "Notch"). You can pass multiple
  pathway names, or `NULL`. By default, all pathways defined in
  [`gnomeR::pathways`](https://mskcc-epi-bio.github.io/gnomeR/reference/pathways.md)
  will be included. Included default pathways are alteration-specific,
  meaning a specific type of alteration (mut/cna/fusion) is required to
  mark a 1 for that pathway.

- custom_pathways:

  a vector of alterations to annotate as a single pathway, or a list of
  custom pathways (see
  [`gnomeR::pathways`](https://mskcc-epi-bio.github.io/gnomeR/reference/pathways.md)
  as example). You must specify the alteration type for each gene using
  `.mut`, `.Amp`, `.Del` suffix, e.g. `c("TP53.mut", "CDKN2A.Amp")`. If
  you wish to count any type of alteration on that gene towards the
  pathway you can use the `.any` suffix (e.g. `c("TP53.any")`).

- other_vars:

  One or more column names (quoted or unquoted) in data to be retained
  in resulting data frame. Default is NULL.

- count_pathways_by:

  deprecated

## Value

a data frame: each sample is a row, columns are pathways, with values of
0/1 depending on pathway alteration status.

## Details

Input a binary matrix of patients x alterations and return a dataframe
with a column per pathway indicating if default or custom oncogenic
signaling pathways are activated in each sample. Default package
pathways were sourced from [Sanchez-Vega, F et al.,
2018](https://pubmed.ncbi.nlm.nih.gov/29625050/).

Please check for gene aliases in your data set before using.

## Examples

``` r
gene_binary <- create_gene_binary(mutation = gnomeR::mutations,
 cna = gnomeR::cna,
 fusion = gnomeR::sv)
#> ! `samples` argument is `NULL`. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in mutation, fusion or cna data frames
#> ! 7 mutations have `NA` or blank in the mutationStatus column instead of 'SOMATIC' or 'GERMLINE'. These were assumed to be 'SOMATIC' and were retained in the resulting binary matrix.
pathway_df <- add_pathways(gene_binary, pathways = "Notch")
```
