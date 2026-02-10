# Check for silent mutations

Check for silent mutations

## Usage

``` r
.check_for_silent(mutation, include_silent)
```

## Arguments

- mutation:

  Raw maf dataframe containing alteration data

- include_silent:

  Silent mutations will be removed if FALSE (default). Variant
  classification column is needed.

## Value

a corrected maf file or an error if problems with maf
