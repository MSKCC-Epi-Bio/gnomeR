# provide a list of impact panels a provided gene is found within

provide a list of impact panels a provided gene is found within

## Usage

``` r
which_impact_panel(hugo_symbol)
```

## Arguments

- hugo_symbol:

  a vector of hugo symbols

## Value

a data frame with hugo symbols and the IMPACT panels on which they are
included

## Examples

``` r
hugos <- unique(gnomeR::mutations$hugoGeneSymbol)[1:10]

which_impact_panel(hugos)
#> # A tibble: 10 × 7
#>    genes_in_panel `IMPACT-HEME-400` `IMPACT-HEME-468` IMPACT341 IMPACT410
#>    <chr>          <chr>             <chr>             <chr>     <chr>    
#>  1 AKT1           yes               yes               yes       yes      
#>  2 ALK            yes               yes               yes       yes      
#>  3 APC            yes               yes               yes       yes      
#>  4 AR             yes               yes               yes       yes      
#>  5 ARAF           yes               yes               yes       yes      
#>  6 ATM            yes               yes               yes       yes      
#>  7 ATR            yes               yes               yes       yes      
#>  8 PARP1          yes               yes               yes       yes      
#>  9 RHOA           yes               yes               yes       yes      
#> 10 ZFHX3          no                no                no        yes      
#> # ℹ 2 more variables: IMPACT468 <chr>, IMPACT505 <chr>
```
