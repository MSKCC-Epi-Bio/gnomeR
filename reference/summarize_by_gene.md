# Simplify binary matrix to one column per gene that counts any alteration type as 1

This will reduce the number of columns in your binary matrix, and the
resulting data frame will have only 1 col per gene, as opposed to
separate columns for mutation/cna/fusion.

## Usage

``` r
summarize_by_gene(gene_binary, other_vars = NULL)
```

## Arguments

- gene_binary:

  a 0/1 matrix of gene alterations

- other_vars:

  One or more column names (quoted or unquoted) in data to be retained
  in resulting data frame. Default is NULL.

## Value

a binary matrix with a row for each sample and one column per gene

## Examples

``` r
samples <- unique(gnomeR::mutations$sampleId)[1:10]
gene_binary <- create_gene_binary(
  samples = samples, mutation = mutations, cna = cna,
  mut_type = "somatic_only",
  include_silent = FALSE,
  specify_panel = "IMPACT341"
) %>%
  summarize_by_gene()
```
