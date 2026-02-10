# Extract IMPACT Patient ID From Sample ID

Extract IMPACT Patient ID From Sample ID

## Usage

``` r
extract_patient_id(sample_id)
```

## Arguments

- sample_id:

  A character vector of IMPACT Tumor sample IDs

## Value

Returns a vector of patient IDs

## Examples

``` r
sample_id = c("P-0000071-T01-IM3", "P-0000072-T02-IM4", "P-0000073-T03-IM5")
extract_patient_id(sample_id)
#> [1] "P-0000071" "P-0000072" "P-0000073"
```
