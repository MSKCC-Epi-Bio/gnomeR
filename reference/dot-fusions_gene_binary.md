# Make Binary Matrix From Fusion data frame

Make Binary Matrix From Fusion data frame

## Usage

``` r
.fusions_gene_binary(fusion, samples, specify_panel)
```

## Arguments

- fusion:

  A data frame of fusions. If not NULL the outcome will be added to the
  matrix with columns ending in ".fus". Default is NULL.

- samples:

  a character vector specifying which samples should be included in the
  resulting data frame. Default is NULL is which case all samples with
  an alteration in the mutation, cna or fusions file will be used. If
  you specify a vector of samples that contain samples not in any of the
  passed genomic data frames, 0's (or NAs when appropriate if specifying
  a panel) will be returned for every column of that patient row.

- specify_panel:

  Default is `"no"` where no panel annotation is done. Otherwise pass a
  character vector of length 1 with a panel id (see
  [`gnomeR::gene_panels`](https://mskcc-epi-bio.github.io/gnomeR/reference/gene_panels.md)
  for available panels), or `"impact"` for automated IMPACT annotation.
  Alternatively, you may pass a data frame of `sample_id`-`panel_id`
  pairs specifying panels for each sample for which to insert NAs
  indicating genes not tested. See below for details.

## Value

a data frame
