# IMPACT Alias Tables

Data frame of genes and their aliases for IMPACT panel genes. This is
used for gene name resolution functionality.

## Usage

``` r
impact_alias_table
```

## Format

A data frame with 1658 rows

- hugo_symbol:

  gene Hugo Symbol

- alias:

  Alias of Hugo Symbol in `hugo_symbol` column

- entrez_id:

  entrez ID of gene in `hugo_symbol`

- alias_entrez_id:

  entrez ID of `alias` gene
