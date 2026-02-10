# Create binary data.frames depending on type of mutation data

Create binary data.frames depending on type of mutation data

## Usage

``` r
.process_binary(data, samples, type = c("mut", "del", "amp", "fus"))
```

## Arguments

- data:

  a dataset of alterations

- samples:

  a vector of unique sample ids

- type:

  a character indicator for which type of alteration the dataset
  contains

## Value

a data.frame of alterations
