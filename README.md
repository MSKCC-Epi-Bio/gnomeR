
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnomeR

<!-- badges: start -->

[![R-CMD-check](https://github.com/MSKCC-Epi-Bio/gnomeR/workflows/R-CMD-check/badge.svg)](https://github.com/MSKCC-Epi-Bio/gnomeR/actions)
[![Codecov test
coverage](https://codecov.io/gh/MSKCC-Epi-Bio/gnomeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/MSKCC-Epi-Bio/gnomeR?branch=main)

<!-- badges: end -->

## Installation

You can install the development version of `gnomeR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MSKCC-Epi-Bio/gnomeR")
```

Along with its companion package for cbioPortal data download:

``` r
devtools::install_github("karissawhiting/cbioportalR")
```

## Introduction

the `gnomeR` package provides a consistent framework for genetic data
processing, visualization and analysis. This is primarily targeted to
IMPACT datasets but can also be applied to any genomic data provided by
cBioPortal. With {gnomeR} and {cbioportalR} you can:

- [**Download and gather data from
  CbioPortal**](https://www.karissawhiting.com/cbioportalR/) - Pull from
  cBioPortal data base by study ID or sample ID.
- **OncoKB annotate data (coming soon)** - Filter genomic data for known
  oncogenic alterations.
- **Process genomic data** - Process retrieved mutation/maf, fusions,
  copy-number alteration, and segmentation data (when available) into an
  analysis-ready formats.
- **Visualize processed data** - Create summary plots from processed
  data.
- **Analyzing processed data**- Analyze associations between genomic
  variables and clinical variables or outcomes.

{gnomeR} is part of
[gnomeverse](https://mskcc-epi-bio.github.io/genomeverse/), a collection
of R packages designed to work together seamlessly to create
reproducible clinico-genomic analysis pipelines.

## Getting Set up

{gnomeR} works with any genomic data that follows cBioPortal guidelines
for
[mutation](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#data-file-5),
[CNA](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#discrete-copy-number-data),
or
[fusion](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#structural-variant-data)
data file formats.

If you wish to pull the data directly from cBioPortal, see how to get
set up with credentials with the
[{cbioportalR}](https://karissawhiting.github.io/cbioportalR/) package.

## Processing Genomic Data

The below examples uses the data sets `mutatations`, `sv`, `cna` which
were pulled from cBioPortal and are included in the package as example
data sets. We will sample 100 samples for examples:

``` r
set.seed(123)

mut <- gnomeR::mutations
cna <- gnomeR::cna
sv <- gnomeR::sv

un <-  unique(mut$sampleId)
sample_patients <- sample(un, size = 50, replace = FALSE)
```

The main data processing function is `create_gene_binary()` which takes
mutation, CNA and fusion files as input, and outputs a binary matrix of
N rows (number of samples) by M genes included in the data set. We can
specify which patients are included which will force all patients in
resulting dataframe, even if they have no alterations.

``` r
gen_dat <- create_gene_binary(samples = sample_patients,
                         mutation = mut,
                         fusion = sv,
                         cna = cna)

head(gen_dat[, 1:6])
#> # A tibble: 6 × 6
#>   sample_id           ALK   APC    AR  ARAF   ATM
#>   <chr>             <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 P-0004508-T01-IM5     1     0     0     0     0
#> 2 P-0005806-T01-IM5     0     1     0     0     0
#> 3 P-0007006-T01-IM5     0     1     0     0     0
#> 4 P-0008682-T01-IM5     0     1     0     0     0
#> 5 P-0001297-T01-IM3     0     0     1     0     0
#> 6 P-0007538-T01-IM5     0     0     0     1     0
```

By default, mutations, CNA and fusions will be returned in separate
columns. You can combine these at the gene level using the following:

``` r
by_gene <- gen_dat %>% 
  summarize_by_gene()

head(by_gene[,1:6])
#> # A tibble: 6 × 6
#>   sample_id           ALK  ARAF   BLM CD79B CSF1R
#>   <chr>             <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 P-0004508-T01-IM5     1     0     0     0     0
#> 2 P-0005806-T01-IM5     0     0     0     0     0
#> 3 P-0007006-T01-IM5     0     0     0     0     0
#> 4 P-0008682-T01-IM5     0     0     0     0     0
#> 5 P-0001297-T01-IM3     0     0     0     0     0
#> 6 P-0007538-T01-IM5     0     1     0     0     1
```

## Visualize

You can visualize your processed and raw alteration data sets using
{gnomeR}’s many data visualization functions.

Quickly visualize mutation characteristics with `ggvarclass()`,
`ggvartype()`, `ggsnvclass()`, `ggsamplevar()`, `ggtopgenes()`,
`gggenecor()`, and `ggcomut()`.

``` r
ggvarclass(mutation = mut)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Summarize & Analyze

You can tabulate summarize your genomic data frame using the
`tbl_genomic()` function, a wrapper for `gtsummary::tbl_summary()`.

``` r

gen_dat <- gen_dat %>%
  dplyr::mutate(trt_status = sample(x = c("pre-trt", "post-trt"),
       size = nrow(gen_dat), replace = TRUE)) 
```

``` r

gene_tbl_trt <-  gen_dat %>%
  subset_by_frequency(t = .1, other_vars = trt_status) %>%
  tbl_genomic(by = trt_status) %>%
  gtsummary::add_p() 
```

<img src="man/figures/README-tbl_genomic_print.png" width="30%" />

Additionally, you can analyze custom pathways, or a set of default gene
pathways using `add_pathways()`:

``` r

path_by_trt <- gen_dat %>%
  add_pathways() %>%
  select(sample_id, trt_status, contains("pathway_")) %>%
  tbl_genomic(by = trt_status) %>%
  gtsummary::add_p() 
```

<img src="man/figures/README-path_by_trt.png" width="30%" />

# Contributing

Please note that the gnomeR project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

Thank you to all contributors!

[@akriti21](https://github.com/akriti21),
[@alrein-05](https://github.com/alrein-05),
[@arorarshi](https://github.com/arorarshi),
[@AxelitoMartin](https://github.com/AxelitoMartin),
[@brombergm](https://github.com/brombergm),
[@carokos](https://github.com/carokos),
[@ChristineZ-msk](https://github.com/ChristineZ-msk),
[@ddsjoberg](https://github.com/ddsjoberg),
[@edrill](https://github.com/edrill),
[@hfuchs5](https://github.com/hfuchs5),
[@jalavery](https://github.com/jalavery),
[@jflynn264](https://github.com/jflynn264),
[@karissawhiting](https://github.com/karissawhiting),
[@michaelcurry1123](https://github.com/michaelcurry1123),
[@mljaniczek](https://github.com/mljaniczek),
[@slb2240](https://github.com/slb2240),
[@stl2137](https://github.com/stl2137),
[@toumban1](https://github.com/toumban1),
[@whitec4](https://github.com/whitec4), and
[@Yukodeng](https://github.com/Yukodeng)

# The End
