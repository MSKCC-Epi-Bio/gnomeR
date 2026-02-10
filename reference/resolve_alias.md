# Resolve Hugo Symbol Names with Aliases

Resolve Hugo Symbol Names with Aliases

## Usage

``` r
resolve_alias(gene_to_check, alias_table)
```

## Arguments

- gene_to_check:

  hugo_symbol to be check

- alias_table:

  table containing all the aliases

## Value

if the accepted hugo symbol is input, it is returned back. If an alias
name is provided, the more common name/more up to date name is returned

## Examples

``` r
resolve_alias("MLL4", alias_table = impact_alias_table)
#> [1] "KMT2B"
```
