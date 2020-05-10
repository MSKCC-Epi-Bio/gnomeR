
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnomeR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/AxelitoMartin/gnomeR.svg?branch=development)](https://travis-ci.org/AxelitoMartin/gnomeR)
[![Codecov test
coverage](https://codecov.io/gh/AxelitoMartin/gnomeR/branch/development/graph/badge.svg)](https://codecov.io/gh/AxelitoMartin/gnomeR?branch=development)
<!-- badges: end -->

the {gnomeR} package provides a consistent framework for genetic data processing, visualization and analysis. This is primarily targeted to IMPACT datasets but can also be applied to any genomic data provided by CbioPortal. 

  - <div class="text-blue mb-2"> Dowloading and gathering data from CbioPortal</div> through an integrated API using simply the sample IDs of the samples of interests or the name of the study to retrive all samples in that study. 
  - <div class="text-blue mb-2"> Processing genomic data</div> retrieved for mutations (MAF file), fusions (MAF file) and copy-number alterations (and when available segmentation files) into an analysis ready format. 
  - <div class="text-blue mb-2"> Visualization of the processed data</div> provided through MAF file summaries, OncoPrints and heatmaps.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AxelitoMartin/gnomeR")
```
