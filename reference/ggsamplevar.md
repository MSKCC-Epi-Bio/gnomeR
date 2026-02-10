# \#' Utility Function to Extract SNV \#' \#' @param x string \#' @param n number of characters from right \#' \#' @return string \#' @noRd \#' @examples \#' substrRight("Hello", 2)

Histogram of Variants Per Sample Colored By Variant Classification

## Usage

``` r
ggsamplevar(mutation)
```

## Arguments

- mutation:

  Raw mutation dataframe containing alteration data

## Value

Histogram of counts of variants per tumor sample

## Examples

``` r
ggsamplevar(gnomeR::mutations)

```
