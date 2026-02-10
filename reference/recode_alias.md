# Recode Hugo Symbol Column

Searches the Hugo Symbol column in a genomic dataframe to look for any
genes that have common gene name aliases, and replaces those aliases
with the accepted (most recent) gene name. Function uses
[`gnomeR::impact_alias_table`](https://mskcc-epi-bio.github.io/gnomeR/reference/impact_alias_table.md)
by default as reference for which aliases to replace and supports IMPACT
panel alias replacement only at this time. Custom tables can be provided
as long as `hugo_symbol` and `alias` columns exist.

## Usage

``` r
recode_alias(genomic_df, alias_table = "impact", supress_warnings = FALSE)
```

## Arguments

- genomic_df:

  a gene_binary object

- alias_table:

  a string indicating "impact" or "genie", or a dataframe with at least
  two columns (`hugo_symbol`, `alias`) with one row for each pair.

- supress_warnings:

  If TRUE, function will return a list containing a dataframe of recoded
  results and a names vector of recoded aliases in data

## Value

A dataframe with recoded Hugo Symbol columns

## Examples

``` r
genomic_df <- rename_columns(gnomeR::mutations[1:5, ])

alias_table <- data.frame("hugo_symbol" = c("New Symbol", "New Symbol2"),
"alias" = c("PARP1", "AKT1"))

recode_alias(genomic_df, alias_table)
#> Warning: To ensure gene with multiple names/aliases are correctly grouped together, the
#> following genes in your dataframe have been recoded (if you are running
#> `create_gene_binary()` you can prevent this with `alias_table = FALSE`):
#> ! PARP1 recoded to New Symbol
#> ! AKT1 recoded to New Symbol2
#> # A tibble: 5 × 29
#>   hugo_symbol entrez_gene_id uniqueSampleKey                    uniquePatientKey
#>   <chr>                <int> <chr>                              <chr>           
#> 1 New Symbol             142 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxMTI4OnB…
#> 2 New Symbol             142 UC0wMDAxODU5LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxODU5OnB…
#> 3 New Symbol             142 UC0wMDAxODk1LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxODk1OnB…
#> 4 New Symbol2            207 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxMTI4OnB…
#> 5 New Symbol2            207 UC0wMDAxODQ1LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxODQ1OnB…
#> # ℹ 25 more variables: molecular_profile_id <chr>, sample_id <chr>,
#> #   patient_id <chr>, study_id <chr>, center <chr>, mutation_status <chr>,
#> #   validation_status <chr>, start_position <int>, end_position <int>,
#> #   reference_allele <chr>, hgv_sp_short <chr>, variant_classification <chr>,
#> #   functionalImpactScore <chr>, fisValue <dbl>, linkXvar <chr>, linkPdb <chr>,
#> #   linkMsa <chr>, ncbi_build <chr>, variant_type <chr>, keyword <chr>,
#> #   chromosome <chr>, allele <chr>, refseqMrnaId <chr>, …
recode_alias(genomic_df, alias_table, supress_warnings = TRUE)
#> $genomic_df
#> # A tibble: 5 × 29
#>   hugo_symbol entrez_gene_id uniqueSampleKey                    uniquePatientKey
#>   <chr>                <int> <chr>                              <chr>           
#> 1 New Symbol             142 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxMTI4OnB…
#> 2 New Symbol             142 UC0wMDAxODU5LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxODU5OnB…
#> 3 New Symbol             142 UC0wMDAxODk1LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxODk1OnB…
#> 4 New Symbol2            207 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxMTI4OnB…
#> 5 New Symbol2            207 UC0wMDAxODQ1LVQwMS1JTTM6cHJhZF9tc… UC0wMDAxODQ1OnB…
#> # ℹ 25 more variables: molecular_profile_id <chr>, sample_id <chr>,
#> #   patient_id <chr>, study_id <chr>, center <chr>, mutation_status <chr>,
#> #   validation_status <chr>, start_position <int>, end_position <int>,
#> #   reference_allele <chr>, hgv_sp_short <chr>, variant_classification <chr>,
#> #   functionalImpactScore <chr>, fisValue <dbl>, linkXvar <chr>, linkPdb <chr>,
#> #   linkMsa <chr>, ncbi_build <chr>, variant_type <chr>, keyword <chr>,
#> #   chromosome <chr>, allele <chr>, refseqMrnaId <chr>, …
#> 
#> $aliases_in_data
#>                             !                             ! 
#> "PARP1 recoded to New Symbol" "AKT1 recoded to New Symbol2" 
#> 
```
