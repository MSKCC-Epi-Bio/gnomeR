# Reformat Wide CNA Data to Long

Takes a numeric vector of CNA data in wide format where each column is a
sample and each row is a hugo symbol. Function will return a long format
CNA data set of just events (neutral/diploid instances are filtered out)
and will recode events from numeric to descriptive (-2/-1/-1.5 is
deletion, 2/1 is amplification).

## Usage

``` r
pivot_cna_longer(wide_cna, clean_sample_ids = TRUE)
```

## Arguments

- wide_cna:

  a cna dataframe in wide format (e.g. gnomeR::cna)

- clean_sample_ids:

  `TRUE` by default and function will clean `sample_id` field to replace
  "." with "-". If `FALSE`, no modification will be made to returned
  `sample_ids` field

## Value

A long data frame of CNA events

## Examples

``` r
cna <- pivot_cna_longer(wide_cna = gnomeR::cna_wide)
#> ! Replacing all `.` to `-` in sample_id field (e.g. `P.0001930.T01.IM3` -> `P-0001930-T01-IM3`).
#> To prevent this, use argument `clean_sample_ids = FALSE`
```
