# Barplot of Most Frequently Altered Genes

Barplot of Most Frequently Altered Genes

## Usage

``` r
ggtopgenes(mutation, n_genes = 10)
```

## Arguments

- mutation:

  Raw mutation dataframe containing alteration data

- n_genes:

  Number of top genes to display in plot

## Value

Barplot of counts of top variant genes

## Examples

``` r
ggtopgenes(gnomeR::mutations)

```
