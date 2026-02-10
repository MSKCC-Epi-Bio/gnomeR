# Enables creation of a binary matrix from a mutation, fusion or CNA file with a predefined list of samples (rows are samples and columns are genes)

Enables creation of a binary matrix from a mutation, fusion or CNA file
with a predefined list of samples (rows are samples and columns are
genes)

## Usage

``` r
create_gene_binary(
  samples = NULL,
  mutation = NULL,
  mut_type = c("omit_germline", "somatic_only", "germline_only", "all"),
  snp_only = FALSE,
  include_silent = FALSE,
  fusion = NULL,
  cna = NULL,
  high_level_cna_only = FALSE,
  specify_panel = "no",
  recode_aliases = "impact"
)
```

## Arguments

- samples:

  a character vector specifying which samples should be included in the
  resulting data frame. Default is NULL is which case all samples with
  an alteration in the mutation, cna or fusions file will be used. If
  you specify a vector of samples that contain samples not in any of the
  passed genomic data frames, 0's (or NAs when appropriate if specifying
  a panel) will be returned for every column of that patient row.

- mutation:

  A data frame of mutations in the format of a maf file.

- mut_type:

  The mutation type to be used. Options are "omit_germline",
  "somatic_only", "germline_only" or "all". Note "all" will keep all
  mutations regardless of status (not recommended). Default is
  omit_germline which includes all somatic mutations, as well as any
  unknown mutation types (most of which are usually somatic)

- snp_only:

  Boolean to rather the genetics events to be kept only to be SNPs
  (insertions and deletions will be removed). Default is FALSE.

- include_silent:

  Boolean to keep or remove all silent mutations. TRUE keeps, FALSE
  removes. Default is FALSE.

- fusion:

  A data frame of fusions. If not NULL the outcome will be added to the
  matrix with columns ending in ".fus". Default is NULL.

- cna:

  A data frame of copy number alterations. If inputed the outcome will
  be added to the matrix with columns ending in ".del" and ".amp".
  Default is NULL.

- high_level_cna_only:

  If TRUE, only deep deletions (-2, -1.5) or high level
  amplifications (2) will be counted as events in the binary matrix.
  Gains (1) and losses (1) will be ignored. Default is `FALSE` where all
  CNA events are counted.

- specify_panel:

  Default is `"no"` where no panel annotation is done. Otherwise pass a
  character vector of length 1 with a panel id (see
  [`gnomeR::gene_panels`](https://mskcc-epi-bio.github.io/gnomeR/reference/gene_panels.md)
  for available panels), or `"impact"` for automated IMPACT annotation.
  Alternatively, you may pass a data frame of `sample_id`-`panel_id`
  pairs specifying panels for each sample for which to insert NAs
  indicating genes not tested. See below for details.

- recode_aliases:

  Default is `"impact"` where function will check for IMPACT genes that
  may go by more than 1 name in your data and replace the alias name
  with the standardized gene name (see
  [`gnomeR::impact_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/impact_alias_table.md)
  for reference list). If `"no"`, no alias annotation will be performed.
  If `"genie"`, an alias table with GENIE BPC genes will be used to
  check (see
  [`gnomeR::genie_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/genie_alias_table.md)
  for reference list). Alternatively, you may pass a custom alias list
  as a data frame with columns `hugo_symbol` and `alias` specifying a
  custom alias table to use for checks. See below for details.

## Value

a data frame with sample_id and alteration binary columns with values of
0/1

## `specify_panel` argument

- If `specify_panel = "no"` is passed (default) data will be returned as
  is without any additional NA annotations.

- If a single panel id is passed (e.g. `specify_panel = "IMPACT468"`),
  all genes in your data that are not tested on that panel will be set
  to `NA` in results for all samples (see gnomeR::gene_panels to see
  which genes are on each supported panels).

- If `specify_panel = "impact"` is passed, impact panel version will be
  inferred based on each sample_id (based on `IMX` nomenclature) and
  NA's will be annotated accordingly for each sample/panel pair.

- If you wish to specify different panels for each sample, pass a data
  frame (with all samples included) with columns: `sample_id`, and
  `panel_id`. Each sample will be annotated with NAs according to that
  specific panel. If a sample in your data is missing from the
  `sample_id` column in the `specify_panel` dataframe, it will be
  returned with no annotation (equivalent of setting it to "no").

## `recode_aliases` argument

- If `recode_aliases = "impact"` is passed (default), function will use
  [`gnomeR::impact_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/impact_alias_table.md)
  to find and replace any non-standard hugo symbol names with their more
  common (or more recent) accepted gene name.

- If `recode_aliases = "genie"` is passed, function will use
  [`gnomeR::genie_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/genie_alias_table.md)
  to find and replace any non-standard hugo symbol names with their more
  common (or more recent) accepted gene name.

- If `recode_aliases = "no"` is passed, data will be returned as is
  without any alias replacements.

- If you have a custom table of vetted aliases you wish to use, you can
  pass a data frame with columns: `hugo_symbol`, and `alias`. Each row
  should have one gene in the `hugo_symbol` column indicating the
  accepted gene name, and one gene in the `alias` column indicating an
  alias you want to check for and replace. If a gene has multiple
  aliases to check for, each should be represented in its own separate
  row. See
  [`gnomeR::impact_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/impact_alias_table.md)
  for an example of accepted data formatting.

## Examples

``` r
mut.only <- create_gene_binary(mutation = gnomeR::mutations)
#> ! `samples` argument is `NULL`. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in mutation, fusion or cna data frames
#> ! 7 mutations have `NA` or blank in the mutationStatus column instead of 'SOMATIC' or 'GERMLINE'. These were assumed to be 'SOMATIC' and were retained in the resulting binary matrix.

samples <- gnomeR::mutations$sampleId

bin.mut <- create_gene_binary(
  samples = samples, mutation = gnomeR::mutations,
  mut_type = "omit_germline", snp_only = FALSE,
  include_silent = FALSE
)
#> ! 7 mutations have `NA` or blank in the mutationStatus column instead of 'SOMATIC' or 'GERMLINE'. These were assumed to be 'SOMATIC' and were retained in the resulting binary matrix.
```
