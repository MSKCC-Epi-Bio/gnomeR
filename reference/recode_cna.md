# Function to recode numeric CNA alteration values to factor values

Function to recode numeric CNA alteration values to factor values

## Usage

``` r
recode_cna(alteration_vector)
```

## Arguments

- alteration_vector:

  a vector of CNA alterations coded with any of the following levels:
  neutral, deletion, amplification, gain, loss, homozygous deletion,
  hemizygous deletion, loh, gain, high level amplification, 0, -1, -1.5,
  -2, 1, 2.

## Value

a recoded CNA data set with factor alteration values. See details for
code dictionary

## Details

CNA is coded to the following key based on key: values below

- "neutral": "0", "neutral",

- "deletion": "homozygous deletion", "-2",

- "deletion": "loh", "-1.5",

- "deletion": "hemizygous deletion", "-1",

- "amplification": "gain", "1",

- "amplification": high level amplification", "2",

## Examples

``` r
recode_cna(gnomeR::cna$alteration[1:10])
#>  [1] deletion      amplification deletion      amplification deletion     
#>  [6] deletion      deletion      deletion      deletion      amplification
#> Levels: deletion amplification
```
