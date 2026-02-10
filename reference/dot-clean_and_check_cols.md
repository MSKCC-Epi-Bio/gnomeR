# Checks genomic input file columns to ensure column names are correct

Checks genomic input file columns to ensure column names are correct

## Usage

``` r
.clean_and_check_cols(
  df_to_check,
  required_cols = c("sample_id", "hugo_symbol")
)
```

## Arguments

- df_to_check:

  Raw maf dataframe containing alteration data

- required_cols:

  A character specifying names of columns to check

## Value

a corrected maf file or an error if problems with maf

## Examples

``` r
gnomeR:::.clean_and_check_cols(df_to_check = gnomeR::mutations)
#> # A tibble: 725 × 29
#>    hugo_symbol entrez_gene_id uniqueSampleKey                   uniquePatientKey
#>    <chr>                <int> <chr>                             <chr>           
#>  1 PARP1                  142 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9t… UC0wMDAxMTI4OnB…
#>  2 PARP1                  142 UC0wMDAxODU5LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODU5OnB…
#>  3 PARP1                  142 UC0wMDAxODk1LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODk1OnB…
#>  4 AKT1                   207 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9t… UC0wMDAxMTI4OnB…
#>  5 AKT1                   207 UC0wMDAxODQ1LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODQ1OnB…
#>  6 AKT1                   207 UC0wMDA1NTcwLVQwMS1JTTU6cHJhZF9t… UC0wMDA1NTcwOnB…
#>  7 ALK                    238 UC0wMDAxNzY4LVQwMS1JTTM6cHJhZF9t… UC0wMDAxNzY4OnB…
#>  8 ALK                    238 UC0wMDA0NTA4LVQwMS1JTTU6cHJhZF9t… UC0wMDA0NTA4OnB…
#>  9 ALK                    238 UC0wMDAxODk1LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODk1OnB…
#> 10 ALK                    238 UC0wMDAyOTg0LVQwMS1JTTM6cHJhZF9t… UC0wMDAyOTg0OnB…
#> # ℹ 715 more rows
#> # ℹ 25 more variables: molecular_profile_id <chr>, sample_id <chr>,
#> #   patient_id <chr>, study_id <chr>, center <chr>, mutation_status <chr>,
#> #   validation_status <chr>, start_position <int>, end_position <int>,
#> #   reference_allele <chr>, hgv_sp_short <chr>, variant_classification <chr>,
#> #   functionalImpactScore <chr>, fisValue <dbl>, linkXvar <chr>, linkPdb <chr>,
#> #   linkMsa <chr>, ncbi_build <chr>, variant_type <chr>, keyword <chr>, …
```
