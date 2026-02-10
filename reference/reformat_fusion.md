# Enables users to reformat fusions files so that each fusion is listed as one row with two hugo-symbol sites instead of two rows, one for each site. This is the required format for the `create_gene_binary` function.

Enables users to reformat fusions files so that each fusion is listed as
one row with two hugo-symbol sites instead of two rows, one for each
site. This is the required format for the `create_gene_binary` function.

## Usage

``` r
reformat_fusion(fusions)
```

## Arguments

- fusions:

  a data frame of fusion/structural variants that occur in a cohort.
  There should be a `sample_id`, `hugo_symbol`, and `fusion` column at
  minimum. Intragenic/intergenic fusions will have one row. Any two gene
  fusions will have two rows. See
  [`gnomeR::sv_long`](https://mskcc-epi-bio.github.io/gnomeR/reference/sv_long.md)
  for an example.

## Value

a data frame with `sample_id`, `site1hugo_symbol`, and
`site2hugo_symbol` and `fusion` columns. This should match the format of
the
[`gnomeR::sv`](https://mskcc-epi-bio.github.io/gnomeR/reference/sv.md)
dataset.

## Examples

``` r
sv_long1 <- gnomeR::sv_long %>%
  rename_columns() %>%
  reformat_fusion()

head(sv_long1)
#> # A tibble: 6 × 6
#>   sample_id             site_1_hugo_symbol site_2_hugo_symbol site_3_hugo_symbol
#>   <chr>                 <chr>              <chr>              <chr>             
#> 1 GENIE-MSK-P-0005931-… MYD88              OXSR1              NA                
#> 2 GENIE-MSK-P-0013935-… RB1                NA                 NA                
#> 3 GENIE-MSK-P-0013835-… ANKRD11            NA                 NA                
#> 4 GENIE-MSK-P-0012525-… CD74               ROS1               NA                
#> 5 GENIE-MSK-P-0014366-… FBN2               RAD50              NA                
#> 6 GENIE-MSK-P-0010247-… CDCA8              FANCA              NA                
#> # ℹ 2 more variables: site_4_hugo_symbol <chr>, fusion <chr>
```
