# Add descriptive endings to hugo symbol names that do not have one already

Add descriptive endings to hugo symbol names that do not have one
already

## Usage

``` r
.paste_endings(names, ending = NULL)
```

## Arguments

- names:

  hugo symbols to check

- ending:

  character ending to add to hugo symbol names without descriptive
  endings. The default is ".mut". If interested in any type of
  alteration, use ".any".

## Value

a vector of hugo symbols where each entry has a descriptive ending from
the following list: ".Amp", ".Del", ".fus", ".cna", ".mut".
