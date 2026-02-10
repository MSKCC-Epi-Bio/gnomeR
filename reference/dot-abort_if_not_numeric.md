# Check if all gene_binary columns except sample_id and other_vars are numeric

Check if all gene_binary columns except sample_id and other_vars are
numeric

## Usage

``` r
.abort_if_not_numeric(alt_data)
```

## Arguments

- alt_data:

  a binary data frame created from
  [`create_gene_binary()`](https://mskcc-epi-bio.github.io/gnomeR/reference/create_gene_binary.md)

## Value

an error message if not all columns are numeric
