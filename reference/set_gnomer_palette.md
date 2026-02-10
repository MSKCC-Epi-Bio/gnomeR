# Set gnomeR color palette

This function sets the gnomeR color palette as the default palette for
all ggplot2 objects. It does so by overriding the following four
functions from the ggplot2 package:
[`scale_color_discrete()`](https://ggplot2.tidyverse.org/reference/scale_colour_discrete.html),
[`scale_fill_discrete()`](https://ggplot2.tidyverse.org/reference/scale_colour_discrete.html),
[`scale_color_continuous()`](https://ggplot2.tidyverse.org/reference/scale_colour_continuous.html),
and
[`scale_fill_continuous()`](https://ggplot2.tidyverse.org/reference/scale_colour_continuous.html),
and places them in the specified environment. A typical workflow would
include this function at the top of a script, and subsequent calls to
[`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html) will
utilize the gnomeR color palette.

## Usage

``` r
set_gnomer_palette(
  palette = c("pancan", "main", "sunset"),
  gradient = c("pancan", "main", "sunset"),
  reverse = FALSE,
  env = rlang::caller_env()
)
```

## Arguments

- palette:

  name of palette in gnomer_palettes, supplied in quotes. Options
  include `"pancan", "main", "sunset"`. Default is `"pancan"`.

- gradient:

  name of gradient palette in `gnomer_palettes`, supplied in quotes.
  Options include `"pancan", "main", "sunset"`. Default is `"pancan"`.

- reverse:

  if set to `TRUE`, will reverse the order of the color palette

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


# setting other gnomeR palettes
set_gnomer_palette(palette = "main", gradient = "sunset")

ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
  geom_point()


ggplot(mtcars, aes(wt, mpg, color = cyl)) +
  geom_point()
#> Warning: The `scale_name` argument of `continuous_scale()` is deprecated as of ggplot2
#> 3.5.0.
```
