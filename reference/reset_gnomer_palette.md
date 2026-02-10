# Reset gnomeR color palette

This function resets the gnomeR color palette back to the ggplot2
default palette for all ggplot2 objects. A typical workflow would
include this after a call to
[`set_gnomer_palette()`](https://mskcc-epi-bio.github.io/gnomeR/reference/set_gnomer_palette.md)
function is no longer needed, and subsequent calls to
[`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html) will
utilize the default color palette from ggplot2.

## Usage

``` r
reset_gnomer_palette(env = rlang::caller_env())
```

## Arguments

- env:

  environment in which palette will take effect. Default is
  [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html).

## Author

Michael Curry

## Examples

``` r
library(ggplot2)

set_gnomer_palette()

ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point()


reset_gnomer_palette()
#default reset
ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point()

```
