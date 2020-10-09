
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnomeR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/AxelitoMartin/gnomeR.svg?branch=development)](https://travis-ci.org/AxelitoMartin/gnomeR)
[![Codecov test
coverage](https://codecov.io/gh/AxelitoMartin/gnomeR/branch/development/graph/badge.svg)](https://codecov.io/gh/AxelitoMartin/gnomeR?branch=development)
<!-- badges: end -->

the `gnomeR` package provides a consistent framework for genetic data
processing, visualization and analysis. This is primarily targeted to
IMPACT datasets but can also be applied to any genomic data provided by
CbioPortal.

  - [**Dowloading and gathering data from
    CbioPortal**](https://github.com/karissawhiting/cbioportalr) through
    an integrated API using simply the sample IDs of the samples of
    interests or the name of the study to retrive all samples in that
    study. A separate package `cbioportalr` was developed independently.
  - [**Processing genomic
    data**](https://axelitomartin.github.io/gnomeR/articles/Data-processing.html)
    retrieved for mutations (MAF file), fusions (MAF file) and
    copy-number alterations (and when available segmentation files) into
    an analysis ready format.
  - [**Visualization of the processed
    data**](https://axelitomartin.github.io/gnomeR/articles/Visualizations.html)
    provided through MAF file summaries, OncoPrints and heatmaps.
  - [**Analyzing the processed
    data**](https://axelitomartin.github.io/gnomeR/articles/Analizing-genomic-data.html)
    for association with binary, continuous and survival outcome.
    Including further visualiztion to improve understanding of the
    results.

## Installation

You can install `gnomeR` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AxelitoMartin/gnomeR")
```

Similarly for those who wish to explore the development version of
`gnomeR`:

``` r
devtools::install_github("AxelitoMartin/gnomeR", ref = "development")
```

Along with its companion package for cbioPortal data download:

``` r
devtools::install_github("karissawhiting/cbioportalr")
```

## Examples

### Setting up the API

In order to download the data from CbioPortal, one must first require a
token from the website [CbioPortal](https://cbioportal.mskcc.org/) wich
will prompt a login page with your MSKCC credentials. Then navigate to
“Web API” in the top bar menu, following this simply download a token
and copy it after running the following command in R:

``` r
usethis::edit_r_environ()
```

And pasting the token you were given in the .Renviron file that was
created and saving after pasting your token.

``` r
CBIOPORTAL_TOKEN = 'YOUR_TOKEN'
```

You can test your connection using:

``` r
get_cbioportal_token()
```

### Retrieving data

Now that the Cbioportal API is set up in your environment, you must
first specify the database of interest (IMPACT or TCGA are the two
available options). Following this one can either sepcify the samples or
study of interest:

``` r
library(gnomeR)
library(cbioportalr)
ids <- as.character(unique(mut$Tumor_Sample_Barcode)[1:100])
df <- get_genetics(sample_ids = ids,database = "msk_impact",
                       mutations = TRUE, fusions = TRUE, cna = TRUE)
```

### Processing the downloaded data

The `binmat()` function is the feature of the data processing of
`gnomeR`. It takes genomic inputs from various sources of CbioPortal
(mutation files, fusion files and copy number raw counts) to give out a
clean binary matrix of n samples by all the events that were found in
the files.

``` r
df.clean <- binmat(maf = df$mut, cna = df$cna)
```

We further included example datasets from the raw dowloaded files on
CbioPortal (`mut`, `fusion`, `cna`) which we will use for the following
examples.

``` r
set.seed(123)
patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 100, replace=FALSE)]

gen.dat <- binmat(patients = patients, maf = mut, fusion = fusion, cna = cna)
gt(gen.dat[1:10,1:10],rownames_to_stub = TRUE)
```

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#uehsjrjlbj .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#uehsjrjlbj .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#uehsjrjlbj .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#uehsjrjlbj .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#uehsjrjlbj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uehsjrjlbj .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#uehsjrjlbj .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#uehsjrjlbj .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#uehsjrjlbj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#uehsjrjlbj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#uehsjrjlbj .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#uehsjrjlbj .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#uehsjrjlbj .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#uehsjrjlbj .gt_from_md > :first-child {
  margin-top: 0;
}

#uehsjrjlbj .gt_from_md > :last-child {
  margin-bottom: 0;
}

#uehsjrjlbj .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#uehsjrjlbj .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#uehsjrjlbj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#uehsjrjlbj .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#uehsjrjlbj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#uehsjrjlbj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#uehsjrjlbj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#uehsjrjlbj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uehsjrjlbj .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#uehsjrjlbj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#uehsjrjlbj .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#uehsjrjlbj .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#uehsjrjlbj .gt_left {
  text-align: left;
}

#uehsjrjlbj .gt_center {
  text-align: center;
}

#uehsjrjlbj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#uehsjrjlbj .gt_font_normal {
  font-weight: normal;
}

#uehsjrjlbj .gt_font_bold {
  font-weight: bold;
}

#uehsjrjlbj .gt_font_italic {
  font-style: italic;
}

#uehsjrjlbj .gt_super {
  font-size: 65%;
}

#uehsjrjlbj .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="uehsjrjlbj" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

FLT4

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

KRAS

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

TP53

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

NF1

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

ARID1A

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

CARD11

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

ARID5B

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

BCOR

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

MLL2

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

SMAD3

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left gt_stub">

P-0010604-T01-IM5

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0002651-T01-IM3

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0000270-T01-IM3

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0002915-T01-IM3

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0011099-T01-IM5

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0000080-T01-IM3

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0001741-T01-IM3

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0003964-T01-IM3

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

1

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0003842-T01-IM5

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

P-0002597-T02-IM5

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

<td class="gt_row gt_right">

0

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

### Visualization

#### MAF

Before we move on to more complex visualizations, we integrate the
`maf_viz()` function to give an overview of the distribution of the
different mutations across the cohort of interest:

``` r
sum.plots <- maf_viz(maf = mut %>% filter(Tumor_Sample_Barcode %in% patients))
sum.plots$topgenes
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
sum.plots$genecomut
```

<img src="man/figures/README-unnamed-chunk-8-2.png" width="100%" />

#### OncoPrints

OncoPrints are a convenient way to display the overall genomic profiles
of samples in the cohort of interest. This is best used for a subset of
genes that are under consideration.

``` r
genes <- c("TP53","PIK3CA","KRAS","TERT","EGFR","FAT","ALK","CDKN2A","CDKN2B")
plot_oncoPrint(gen.dat = gen.dat %>% select(starts_with(genes)))
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

#### FACETs

[FACETs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5027494/) is an
ASCN tool and open-source software with a broad application to whole
genome, whole-exome, as well as targeted panel sequencing platforms. It
is a fully integrated stand-alone pipeline that includes sequencing BAM
file post-processing, joint segmentation of total- and allele-specific
read counts, and integer copy number calls corrected for tumor purity,
ploidy and clonal heterogeneity, with comprehensive output.

``` r
p.heat <- facets.heatmap(seg = seg, patients = patients, min.purity = 0)
p.heat$p
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

### Analysis

In this section we will quickly overview the possible analysis in
gnomeR.

#### Binary and continuous outcomes

The `gen.summary()` function let’s the user perform a large scale
association between the genomic features present in the `binmat()`
function output and an outcome of choice:

  - binary (unpaired test using Fisher’s exact test and paired test
    using McNemmar’s exact test)
  - continuous (using simple linear regression)

<!-- end list -->

``` r
outcome <- factor(rbinom(n = length(patients),size = 1,prob = 1/2),levels = c("0","1"))
out <- gen.summary(gen.dat = gen.dat,outcome = outcome,filter = 0.05)
gt(out$fits[1:10,],rownames_to_stub = TRUE)
```

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#jxvwzkaysj .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#jxvwzkaysj .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jxvwzkaysj .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#jxvwzkaysj .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#jxvwzkaysj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jxvwzkaysj .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jxvwzkaysj .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#jxvwzkaysj .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#jxvwzkaysj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#jxvwzkaysj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#jxvwzkaysj .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#jxvwzkaysj .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#jxvwzkaysj .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#jxvwzkaysj .gt_from_md > :first-child {
  margin-top: 0;
}

#jxvwzkaysj .gt_from_md > :last-child {
  margin-bottom: 0;
}

#jxvwzkaysj .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#jxvwzkaysj .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#jxvwzkaysj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jxvwzkaysj .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#jxvwzkaysj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jxvwzkaysj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#jxvwzkaysj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#jxvwzkaysj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jxvwzkaysj .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jxvwzkaysj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#jxvwzkaysj .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jxvwzkaysj .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#jxvwzkaysj .gt_left {
  text-align: left;
}

#jxvwzkaysj .gt_center {
  text-align: center;
}

#jxvwzkaysj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#jxvwzkaysj .gt_font_normal {
  font-weight: normal;
}

#jxvwzkaysj .gt_font_bold {
  font-weight: bold;
}

#jxvwzkaysj .gt_font_italic {
  font-style: italic;
}

#jxvwzkaysj .gt_super {
  font-size: 65%;
}

#jxvwzkaysj .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="jxvwzkaysj" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Feature

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Overall

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

0(N=48)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

1(N=52)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

OddsRatio

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Pvalue

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

FDR

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Lower

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Upper

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left gt_stub">

PIK3CA

</td>

<td class="gt_row gt_left">

PIK3CA

</td>

<td class="gt_row gt_left">

12%

</td>

<td class="gt_row gt_left">

18.75%

</td>

<td class="gt_row gt_left">

5.77%

</td>

<td class="gt_row gt_left">

0.27

</td>

<td class="gt_row gt_left">

6.47e-02

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0.04

</td>

<td class="gt_row gt_left">

1.17

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

CDH1

</td>

<td class="gt_row gt_left">

CDH1

</td>

<td class="gt_row gt_left">

6%

</td>

<td class="gt_row gt_left">

10.42%

</td>

<td class="gt_row gt_left">

1.92%

</td>

<td class="gt_row gt_left">

0.17

</td>

<td class="gt_row gt_left">

1.02e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0

</td>

<td class="gt_row gt_left">

1.61

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

PTPRD

</td>

<td class="gt_row gt_left">

PTPRD

</td>

<td class="gt_row gt_left">

8%

</td>

<td class="gt_row gt_left">

12.5%

</td>

<td class="gt_row gt_left">

3.85%

</td>

<td class="gt_row gt_left">

0.28

</td>

<td class="gt_row gt_left">

1.49e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0.03

</td>

<td class="gt_row gt_left">

1.69

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

TERT

</td>

<td class="gt_row gt_left">

TERT

</td>

<td class="gt_row gt_left">

15%

</td>

<td class="gt_row gt_left">

20.83%

</td>

<td class="gt_row gt_left">

9.62%

</td>

<td class="gt_row gt_left">

0.41

</td>

<td class="gt_row gt_left">

1.62e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0.1

</td>

<td class="gt_row gt_left">

1.45

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

EPHA5

</td>

<td class="gt_row gt_left">

EPHA5

</td>

<td class="gt_row gt_left">

5%

</td>

<td class="gt_row gt_left">

8.33%

</td>

<td class="gt_row gt_left">

1.92%

</td>

<td class="gt_row gt_left">

0.22

</td>

<td class="gt_row gt_left">

1.92e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0

</td>

<td class="gt_row gt_left">

2.32

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

TSC1

</td>

<td class="gt_row gt_left">

TSC1

</td>

<td class="gt_row gt_left">

5%

</td>

<td class="gt_row gt_left">

8.33%

</td>

<td class="gt_row gt_left">

1.92%

</td>

<td class="gt_row gt_left">

0.22

</td>

<td class="gt_row gt_left">

1.92e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0

</td>

<td class="gt_row gt_left">

2.32

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

NOTCH3

</td>

<td class="gt_row gt_left">

NOTCH3

</td>

<td class="gt_row gt_left">

5%

</td>

<td class="gt_row gt_left">

8.33%

</td>

<td class="gt_row gt_left">

1.92%

</td>

<td class="gt_row gt_left">

0.22

</td>

<td class="gt_row gt_left">

1.92e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0

</td>

<td class="gt_row gt_left">

2.32

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

CDKN2A.Del

</td>

<td class="gt_row gt_left">

CDKN2A.Del

</td>

<td class="gt_row gt_left">

8%

</td>

<td class="gt_row gt_left">

4.17%

</td>

<td class="gt_row gt_left">

11.54%

</td>

<td class="gt_row gt_left">

2.97

</td>

<td class="gt_row gt_left">

2.72e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0.5

</td>

<td class="gt_row gt_left">

31.58

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

EGFR

</td>

<td class="gt_row gt_left">

EGFR

</td>

<td class="gt_row gt_left">

9%

</td>

<td class="gt_row gt_left">

12.5%

</td>

<td class="gt_row gt_left">

5.77%

</td>

<td class="gt_row gt_left">

0.43

</td>

<td class="gt_row gt_left">

3.05e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0.07

</td>

<td class="gt_row gt_left">

2.17

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

MLL2

</td>

<td class="gt_row gt_left">

MLL2

</td>

<td class="gt_row gt_left">

5%

</td>

<td class="gt_row gt_left">

2.08%

</td>

<td class="gt_row gt_left">

7.69%

</td>

<td class="gt_row gt_left">

3.87

</td>

<td class="gt_row gt_left">

3.64e-01

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

<td class="gt_row gt_left">

0.37

</td>

<td class="gt_row gt_left">

196.71

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

``` r
out$forest.plot
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

``` r
out$vPlot
```

<img src="man/figures/README-unnamed-chunk-11-2.png" width="100%" />

#### Survival analysis

Similarly we include simple tools to perform univariate Cox’s
proportional regression adjusted for false discovery rate in the
`uni.cox()` function.

``` r
time <- rexp(length(patients))
status <- outcome
surv.dat <- as.data.frame(cbind(time,status))
out <- uni.cox(X = gen.dat, surv.dat = surv.dat, surv.formula = Surv(time,status)~.,filter = 0.05)
kable(out$tab[1:10,],row.names = T)
```

|    | Feature    | Coefficient |    HR | Pvalue |       FDR | MutationFrequency |
| :- | :--------- | ----------: | ----: | -----: | --------: | ----------------: |
| 1  | STK11      |       0.930 | 2.540 | 0.0516 | 0.7208300 |              0.07 |
| 2  | CDKN2A.Del |       0.816 | 2.260 | 0.0637 | 0.7208300 |              0.08 |
| 3  | TERT       |     \-0.862 | 0.422 | 0.0703 | 0.7208300 |              0.15 |
| 4  | MLL2       |       0.915 | 2.500 | 0.0859 | 0.7208300 |              0.05 |
| 5  | PTPRD      |     \-1.180 | 0.309 | 0.1040 | 0.7208300 |              0.08 |
| 6  | TSC1       |     \-1.610 | 0.199 | 0.1110 | 0.7208300 |              0.05 |
| 7  | CDH1       |     \-1.430 | 0.238 | 0.1560 | 0.8462973 |              0.06 |
| 8  | PIK3CA     |     \-0.784 | 0.457 | 0.1890 | 0.8462973 |              0.12 |
| 9  | EPHA5      |     \-1.230 | 0.292 | 0.2240 | 0.8462973 |              0.05 |
| 10 | FLT4       |     \-0.870 | 0.419 | 0.2300 | 0.8462973 |              0.05 |

``` r
# out$p
out$KM[[1]]
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

### Further analytical tools

The primary goal of `gnomeR` not being in depth analysis of genomic data
but rather reliable, modulable and reproducible framework for processing
various types of genomic data. For users interested in large scale
genomic analytical methods we compiled various packages developed by
Memorial Sloan-Kettering Cancer Center employees under an umbrella R
package, [gnomeVerse](https://github.com/AxelitoMartin/genomeVerse).
