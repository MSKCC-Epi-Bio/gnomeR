# Subset a Binary Matrix By Alteration Frequency Threshold

Subset a Binary Matrix By Alteration Frequency Threshold

## Usage

``` r
subset_by_frequency(gene_binary, t = 0.1, other_vars = NULL, by = NULL)
```

## Arguments

- gene_binary:

  A data frame with a row for each sample and column for each
  alteration. Data frame must have a `sample_id` column and columns for
  each alteration with values of 0, 1 or NA.

- t:

  Threshold value between 0 and 1 to subset by. Default is 10% (.1).

- other_vars:

  One or more column names (quoted or unquoted) in data to be retained
  in resulting data frame. Default is NULL.

- by:

  Variable used to subset the data. Default is NULL.

## Value

a data frame with a `sample_id` column and columns for alterations over
the given prevalence threshold of `t`.

## Examples

``` r
samples <- unique(gnomeR::mutations$sampleId)
 gene_binary <- create_gene_binary(
   samples = samples, mutation = mutations, cna = cna,
   mut_type = "somatic_only",
   include_silent = FALSE,
   specify_panel = "impact"
 )
gene_binary %>%
 subset_by_frequency()
#> # A tibble: 200 × 7
#>    sample_id          TP53 FOXA1 AR.Amp  SPOP PTEN.Del KMT2C
#>    <chr>             <dbl> <dbl>  <dbl> <dbl>    <dbl> <dbl>
#>  1 P-0001128-T01-IM3     1     0      0     0        0     0
#>  2 P-0001859-T01-IM3     0     0      1     0        0     0
#>  3 P-0001895-T01-IM3     0     0      0     1        0     0
#>  4 P-0001845-T01-IM3     0     0      0     0        0     1
#>  5 P-0001768-T01-IM3     1     0      0     0        0     0
#>  6 P-0002984-T01-IM3     1     0      1     0        0     0
#>  7 P-0000964-T02-IM3     0     0      0     0        0     0
#>  8 P-0000964-T01-IM3     0     0      0     0        0     0
#>  9 P-0000610-T01-IM3     1     0      0     0        1     0
#> 10 P-0001247-T01-IM3     1     0      1     0        0     0
#> # ℹ 190 more rows
```
