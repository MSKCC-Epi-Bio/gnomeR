# Complete list of gnomeR color palettes

Creates color palettes based on gnomeR colors, including "main" which is
a collection of 33 pale-yet-distinct colors which flow nicely (good for
clustering solutions, for example), and "pancan" which is a collection
of 33 vibrant and distinct colors, good for visualizing the 33 TCGA
PanCan cancer types.

## Usage

``` r
gnomer_palettes
```

## Format

An object of class `list` of length 3.

## Examples

``` r
gnomer_palettes[["pancan"]]
#>  [1] "#67000D" "#A50F15" "#EF3B2C" "#FC9272" "#FEE0D2" "#BCBDDC" "#807DBA"
#>  [8] "#54278F" "#3F007D" "#08306B" "#08519C" "#4292C6" "#9ECAE1" "#000000"
#> [15] "#525252" "#969696" "#BDBDBD" "#D9D9D9" "#80CDC1" "#35978F" "#01665E"
#> [22] "#006D2C" "#41AB5D" "#A1D99B" "#FFFFCC" "#FED976" "#FD8D3C" "#8C510A"
#> [29] "#BF812D" "#DFC27D" "#FA9FB5" "#F768A1" "#DD3497"
```
