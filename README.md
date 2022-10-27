
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnomeR

<!-- badges: start -->

[![R-CMD-check](https://github.com/MSKCC-Epi-Bio/gnomeR/workflows/R-CMD-check/badge.svg)](https://github.com/MSKCC-Epi-Bio/gnomeR/actions)
[![Codecov test
coverage](https://codecov.io/gh/MSKCC-Epi-Bio/gnomeR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/MSKCC-Epi-Bio/gnomeR?branch=main)

<!-- badges: end -->

<font size="5">:bangbang: :warning: **NOTE: This package is currently
under active development with a new stable release expected April 30th,
2022. For code written before 2022-03-23, please use the previous stable
version (v1.1.0)**:warning::bangbang: </font>

You can install the pre-2022-03-23 version with:

``` r
remotes::install_github('MSKCC-Epi-Bio/gnomeR@v1.1.0')
```

## Installation

You can install the development version of `gnomeR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MSKCC-Epi-Bio/gnomeR")
```

Along with its companion package for cbioPortal data download:

``` r
devtools::install_github("karissawhiting/cbioportalr")
```

## Introduction

the `gnomeR` package provides a consistent framework for genetic data
processing, visualization and analysis. This is primarily targeted to
IMPACT datasets but can also be applied to any genomic data provided by
cBioPortal. With {gnomeR} and {cbioportalR} you can:

- [**Download and gather data from
  CbioPortal**](https://www.karissawhiting.com/cbioportalR/) - Pull from
  cBioPortal data base by study ID or sample ID.
- **OncoKB annotate data** - Filter genomic data for known oncogenic
  alterations.
- **Process genomic data** - Process retrieved mutation/maf, fusions,
  copy-number alteration, and segmentation data (when available) into an
  analysis-ready formats.
- **Visualize processed data** - Create OncoPrints, heatmaps and summary
  plots from processed data.
- **Analyzing processed data**- Analyze associations between genomic
  variables and clinical variables or outcomes with summary tables,
  advanced visualizations, and univariate models.

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
sample_patients <- sample(un, size = 100, replace = FALSE)
```

The main data processing function is `create_gene_binary()` which takes
mutation, CNA and fusion files as input, and outputs a binary matrix of
N rows (number of samples) by M genes included in the data set. We can
specify which patients are included which will force all patients in
resulting dataframe, even if they have no alterations.

``` r
gen_dat <- create_gene_binary(samples = sample_patients,
                         maf = mut,
                         fusion = sv,
                         cna = cna)

head(gen_dat[, 1:6])
#>                   ERG.fus KDM5A.fus KDM5D.fus GSK3B.fus EGFR.fus PBRM1.fus
#> P-0008869-T01-IM5       1         0         0         0        0         0
#> P-0001242-T01-IM3       0         0         0         0        0         0
#> P-0005806-T01-IM5       0         0         0         0        0         0
#> P-0007346-T01-IM5       0         0         0         0        0         0
#> P-0001861-T01-IM3       0         0         0         0        0         0
#> P-0001202-T01-IM3       0         0         0         0        0         0
```

By default, mutations, CNA and fusions will be returned in separate
columns. You can combine these at teh gene level using the following:

``` r
by_gene <- gen_dat %>% 
  summarize_by_gene()

head(by_gene[,1:6])
#> # A tibble: 6 × 6
#>   sample_id           ERG KDM5A KDM5D GSK3B  EGFR
#>   <chr>             <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 P-0008869-T01-IM5     1     0     0     0     0
#> 2 P-0001242-T01-IM3     0     0     0     0     0
#> 3 P-0005806-T01-IM5     0     0     0     0     0
#> 4 P-0007346-T01-IM5     0     0     0     0     0
#> 5 P-0001861-T01-IM3     0     0     0     0     0
#> 6 P-0001202-T01-IM3     0     0     0     0     0
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

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

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
  tbl_genomic(freq_cutoff = .1, by = trt_status) %>%
  gtsummary::add_p() 

gene_tbl_trt
```

<div id="zdmnupwahh" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#zdmnupwahh .gt_table {
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

#zdmnupwahh .gt_heading {
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

#zdmnupwahh .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#zdmnupwahh .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#zdmnupwahh .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zdmnupwahh .gt_col_headings {
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

#zdmnupwahh .gt_col_heading {
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

#zdmnupwahh .gt_column_spanner_outer {
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

#zdmnupwahh .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#zdmnupwahh .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#zdmnupwahh .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#zdmnupwahh .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
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

#zdmnupwahh .gt_empty_group_heading {
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

#zdmnupwahh .gt_from_md > :first-child {
  margin-top: 0;
}

#zdmnupwahh .gt_from_md > :last-child {
  margin-bottom: 0;
}

#zdmnupwahh .gt_row {
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

#zdmnupwahh .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#zdmnupwahh .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#zdmnupwahh .gt_row_group_first td {
  border-top-width: 2px;
}

#zdmnupwahh .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#zdmnupwahh .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#zdmnupwahh .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#zdmnupwahh .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zdmnupwahh .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#zdmnupwahh .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#zdmnupwahh .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#zdmnupwahh .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#zdmnupwahh .gt_footnotes {
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

#zdmnupwahh .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#zdmnupwahh .gt_sourcenotes {
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

#zdmnupwahh .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#zdmnupwahh .gt_left {
  text-align: left;
}

#zdmnupwahh .gt_center {
  text-align: center;
}

#zdmnupwahh .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#zdmnupwahh .gt_font_normal {
  font-weight: normal;
}

#zdmnupwahh .gt_font_bold {
  font-weight: bold;
}

#zdmnupwahh .gt_font_italic {
  font-style: italic;
}

#zdmnupwahh .gt_super {
  font-size: 65%;
}

#zdmnupwahh .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#zdmnupwahh .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#zdmnupwahh .gt_indent_1 {
  text-indent: 5px;
}

#zdmnupwahh .gt_indent_2 {
  text-indent: 10px;
}

#zdmnupwahh .gt_indent_3 {
  text-indent: 15px;
}

#zdmnupwahh .gt_indent_4 {
  text-indent: 20px;
}

#zdmnupwahh .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>Overall</strong>, N = 100<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>post-trt</strong>, N = 56<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>pre-trt</strong>, N = 44<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col"><strong>p-value</strong><sup class="gt_footnote_marks">2</sup></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left" style="font-weight: bold;">ERG.Amp</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">AR.Amp</td>
<td class="gt_row gt_center">14 (14%)</td>
<td class="gt_row gt_center">6 (11%)</td>
<td class="gt_row gt_center">8 (18%)</td>
<td class="gt_row gt_center">0.3</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">AR.Del</td>
<td class="gt_row gt_center">1 (1.0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (2.3%)</td>
<td class="gt_row gt_center">0.4</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">PTEN.Del</td>
<td class="gt_row gt_center">12 (12%)</td>
<td class="gt_row gt_center">4 (7.1%)</td>
<td class="gt_row gt_center">8 (18%)</td>
<td class="gt_row gt_center">0.092</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">ERG.fus</td>
<td class="gt_row gt_center">27 (27%)</td>
<td class="gt_row gt_center">13 (23%)</td>
<td class="gt_row gt_center">14 (32%)</td>
<td class="gt_row gt_center">0.3</td></tr>
    <tr><td class="gt_row gt_left" style="font-weight: bold;">PTEN.fus</td>
<td class="gt_row gt_center">1 (1.0%)</td>
<td class="gt_row gt_center">0 (0%)</td>
<td class="gt_row gt_center">1 (2.3%)</td>
<td class="gt_row gt_center">0.4</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks">1</sup> n (%)</td>
    </tr>
    <tr>
      <td class="gt_footnote" colspan="5"><sup class="gt_footnote_marks">2</sup> Pearson's Chi-squared test; Fisher's exact test</td>
    </tr>
  </tfoot>
</table>
</div>

Additionally, you can analyze custom pathways, or a set of default gene
pathways using `add_pathways()`:

``` r
path_by_trt <- gen_dat %>%
  add_pathways() %>%
  select(trt_status, contains("pathway_")) %>%
  tbl_genomic(by = trt_status) %>%
  gtsummary::add_p() 


path_by_trt
```

## Further analytical tools

Along with mutation, cna and fusion data, {gnomeR} also allows analysis
and visualization of
[FACETs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5027494/) data.
[FACETs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5027494/) is an
allele-specific copy number tool and open-source software with a broad
application to whole genome, whole-exome, as well as targeted panel
sequencing platforms. It is a fully integrated stand-alone pipeline that
includes sequencing BAM file post-processing, joint segmentation of
total- and allele-specific read counts, and integer copy number calls
corrected for tumor purity, ploidy and clonal heterogeneity, with
comprehensive output.

You can visualize this data using `facets_heatmap()`

``` r

select_samples <- sample(unique(seg$ID),  100)
p.heat <- facets_heatmap(seg = seg, samples = select_samples, min_purity = 0)
p.heat$p
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

# Contributing

Please note that the gnomeR project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

Thank you to all contributors!

[@akriti21](https://github.com/akriti21),
[@arorarshi](https://github.com/arorarshi),
[@AxelitoMartin](https://github.com/AxelitoMartin),
[@carokos](https://github.com/carokos),
[@ChristineZ-msk](https://github.com/ChristineZ-msk),
[@edrill](https://github.com/edrill),
[@jalavery](https://github.com/jalavery),
[@jflynn264](https://github.com/jflynn264),
[@karissawhiting](https://github.com/karissawhiting),
[@michaelcurry1123](https://github.com/michaelcurry1123),
[@mljaniczek](https://github.com/mljaniczek), and
[@slb2240](https://github.com/slb2240)
