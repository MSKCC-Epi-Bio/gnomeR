# Comutation Heatmap of the Top Altered Genes

Comutation Heatmap of the Top Altered Genes

## Usage

``` r
ggcomut(mutation, n_genes = 10, ...)
```

## Arguments

- mutation:

  Raw mutation dataframe containing alteration data

- n_genes:

  Number of top genes to display in plot

- ...:

  Further create_gene_binary() arguments

## Value

Comutation heatmap of the top genes

## Examples

``` r
ggcomut(mutation = gnomeR::mutations)
#> ! `samples` argument is `NULL`. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in mutation, fusion or cna data frames
#> ! 7 mutations have `NA` or blank in the mutationStatus column instead of 'SOMATIC' or 'GERMLINE'. These were assumed to be 'SOMATIC' and were retained in the resulting binary matrix.

```
