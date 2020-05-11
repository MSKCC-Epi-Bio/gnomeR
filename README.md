
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gnomeR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/AxelitoMartin/gnomeR.svg?branch=development)](https://travis-ci.org/AxelitoMartin/gnomeR)
[![Codecov test
coverage](https://codecov.io/gh/AxelitoMartin/gnomeR/branch/development/graph/badge.svg)](https://codecov.io/gh/AxelitoMartin/gnomeR?branch=development)
<!-- badges: end -->

the {gnomeR} package provides a consistent framework for genetic data
processing, visualization and analysis. This is primarily targeted to
IMPACT datasets but can also be applied to any genomic data provided by
CbioPortal.

  - 
    
    <div class="text-blue mb-2">
    
    Dowloading and gathering data from CbioPortal
    
    </div>
    
    through an integrated API using simply the sample IDs of the samples
    of interests or the name of the study to retrive all samples in that
    study.

  - 
    
    <div class="text-blue mb-2">
    
    Processing genomic data
    
    </div>
    
    retrieved for mutations (MAF file), fusions (MAF file) and
    copy-number alterations (and when available segmentation files) into
    an analysis ready format.

  - 
    
    <div class="text-blue mb-2">
    
    Visualization of the processed data
    
    </div>
    
    provided through MAF file summaries, OncoPrints and heatmaps.

  - 
    
    <div class="text-blue mb-2">
    
    Analyzing the processed data
    
    </div>
    
    for association with binary, continuous and survival outcome.
    Including further visualiztion to improve understanding of the
    results.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AxelitoMartin/gnomeR")
```

## Examples

### Setting up the API

In order to download the data from CbioPortal, one must first require a
token from the website [CbioPortal](https://cbiologin.mskcc.org/) wich
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
mycgds = CGDS("https://cbioportal.mskcc.org/",
              token = YOUR_TOKEN)
test(mycgds)
```

### Retrieving data

Now that the Cbioportal API is set up in your environment, you must
first specify the database of interest (IMPACT or TCGA are the two
available options). Following this one can either sepcify the samples or
study of interest:

``` r
get_cbioportal_db("msk_impact")
ids <- as.character(unique(mut$Tumor_Sample_Barcode)[1:100])
df <- get_mutations(sample_ids = ids)

## match genes ##
# This will be included in the future #
df.genes <- df %>% 
  mutate(Hugo_Symbol = as.character(info_impact$Hugo_Symbol[match(entrezGeneId,
                                                                  info_impact$Entrez_Gene_Id)]))
```

### Processing the downloaded data

The `binmat()` function is the feature of the data processing of
`gnomeR`. It takes genomic inputs from various sources of CbioPortal
(mutation files, fusion files and copy number raw counts) to give out a
clean binary matrix of n samples by all the events that were found in
the files.

Still need to deal with stuff for the API here (CNA mostly):

``` r
test2 <- df.genes %>% 
  binmat(maf = ., col.names = c(Tumor_Sample_Barcode = "sampleId", Hugo_Symbol = NULL, Variant_Classification
                                = "mutationType", Mutation_Status = "mutationStatus", Variant_Type = "variantType"))
```

Example with example data in the package for now:

``` r
set.seed(123)
patients <- as.character(unique(mut$Tumor_Sample_Barcode))[sample(1:length(unique(mut$Tumor_Sample_Barcode)), 100, replace=FALSE)]

gen.dat <- binmat(patients = patients, maf = mut, fusion = fusion, cna = cna)
gt(gen.dat[1:10,1:10],rownames_to_stub = T)
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

#uehsjrjlbj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
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
`maf.summary()` function to give an overview of the distribution of the
different mutations across the cohort of interest:

``` r
sum.plots <- maf.summary(maf = mut %>% filter(Tumor_Sample_Barcode %in% patients))
sum.plots$p.genes
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
sum.plots$p.comut
```

<img src="man/figures/README-unnamed-chunk-8-2.png" width="100%" />

#### OncoPrints

OncoPrints are a convenient way to display the overall genomic profiles
of samples in the cohort of interest. This is best used for a subset of
genes that are under consideration.

``` r
genes <- c("TP53","PIK3CA","KRAS","TERT","EGFR","FAT","ALK","CDKN2A","CDKN2B")
plot_oncoPrint(gen.dat = gen.dat %>% select(starts_with(genes)))
#> All mutation types: MUT, AMP, DEL, FUS
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
# need to find how to get facets from API #
# also makes me think we need to get the clinical data too in that for purity #
p.heat <- facets.heatmap(seg = seg, patients = patients, min.purity = 0)
#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion

#> Warning: NAs introduced by coercion
#> The "ward" method has been renamed to "ward.D"; note new "ward.D2"
p.heat$p
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

``` r
# Also need a better example lol#
```

### Analysis

In this section we will quickly overview the possible analysis in
gnomeR.

#### Binary and continuous outcomes

The `gen.tab()` function let’s the user perform a large scale
association between the genomic features present in the `binat()`
function output and an outcome of choice:

  - binary (unpaired test using Fisher’s exact test and paired test
    using McNemmar’s exact test)
  - continuous (using simple linear regression)

<!-- end list -->

``` r
outcome <- factor(rbinom(n = length(patients),size = 1,prob = 1/2),levels = c("0","1"))
out <- gen.tab(gen.dat = gen.dat,outcome = outcome,filter = 0.05)
#> Warning: geom_hline(): Ignoring `mapping` because `yintercept` was provided.
gt(out$fits[1:10,],rownames_to_stub = T)
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

#jxvwzkaysj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
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

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

Overall

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

0(N=48)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

1(N=52)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

OddsRatio

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

Pvalue

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

Lower

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1">

Upper

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

FDR

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left gt_stub">

PIK3CA

</td>

<td class="gt_row gt_center">

12%

</td>

<td class="gt_row gt_center">

18.75%

</td>

<td class="gt_row gt_center">

5.77%

</td>

<td class="gt_row gt_center">

0.27

</td>

<td class="gt_row gt_center">

6.47e-02

</td>

<td class="gt_row gt_center">

0.04

</td>

<td class="gt_row gt_center">

1.17

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

CDH1

</td>

<td class="gt_row gt_center">

6%

</td>

<td class="gt_row gt_center">

10.42%

</td>

<td class="gt_row gt_center">

1.92%

</td>

<td class="gt_row gt_center">

0.17

</td>

<td class="gt_row gt_center">

1.02e-01

</td>

<td class="gt_row gt_center">

0

</td>

<td class="gt_row gt_center">

1.61

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

PTPRD

</td>

<td class="gt_row gt_center">

8%

</td>

<td class="gt_row gt_center">

12.5%

</td>

<td class="gt_row gt_center">

3.85%

</td>

<td class="gt_row gt_center">

0.28

</td>

<td class="gt_row gt_center">

1.49e-01

</td>

<td class="gt_row gt_center">

0.03

</td>

<td class="gt_row gt_center">

1.69

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

TERT

</td>

<td class="gt_row gt_center">

15%

</td>

<td class="gt_row gt_center">

20.83%

</td>

<td class="gt_row gt_center">

9.62%

</td>

<td class="gt_row gt_center">

0.41

</td>

<td class="gt_row gt_center">

1.62e-01

</td>

<td class="gt_row gt_center">

0.1

</td>

<td class="gt_row gt_center">

1.45

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

EPHA5

</td>

<td class="gt_row gt_center">

5%

</td>

<td class="gt_row gt_center">

8.33%

</td>

<td class="gt_row gt_center">

1.92%

</td>

<td class="gt_row gt_center">

0.22

</td>

<td class="gt_row gt_center">

1.92e-01

</td>

<td class="gt_row gt_center">

0

</td>

<td class="gt_row gt_center">

2.32

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

TSC1

</td>

<td class="gt_row gt_center">

5%

</td>

<td class="gt_row gt_center">

8.33%

</td>

<td class="gt_row gt_center">

1.92%

</td>

<td class="gt_row gt_center">

0.22

</td>

<td class="gt_row gt_center">

1.92e-01

</td>

<td class="gt_row gt_center">

0

</td>

<td class="gt_row gt_center">

2.32

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

NOTCH3

</td>

<td class="gt_row gt_center">

5%

</td>

<td class="gt_row gt_center">

8.33%

</td>

<td class="gt_row gt_center">

1.92%

</td>

<td class="gt_row gt_center">

0.22

</td>

<td class="gt_row gt_center">

1.92e-01

</td>

<td class="gt_row gt_center">

0

</td>

<td class="gt_row gt_center">

2.32

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

MYC.Amp

</td>

<td class="gt_row gt_center">

6%

</td>

<td class="gt_row gt_center">

2.08%

</td>

<td class="gt_row gt_center">

9.62%

</td>

<td class="gt_row gt_center">

4.93

</td>

<td class="gt_row gt_center">

2.07e-01

</td>

<td class="gt_row gt_center">

0.52

</td>

<td class="gt_row gt_center">

241.13

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

EGFR

</td>

<td class="gt_row gt_center">

9%

</td>

<td class="gt_row gt_center">

12.5%

</td>

<td class="gt_row gt_center">

5.77%

</td>

<td class="gt_row gt_center">

0.43

</td>

<td class="gt_row gt_center">

3.05e-01

</td>

<td class="gt_row gt_center">

0.07

</td>

<td class="gt_row gt_center">

2.17

</td>

<td class="gt_row gt_left">

9.51e-01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

MLL2

</td>

<td class="gt_row gt_center">

5%

</td>

<td class="gt_row gt_center">

2.08%

</td>

<td class="gt_row gt_center">

7.69%

</td>

<td class="gt_row gt_center">

3.87

</td>

<td class="gt_row gt_center">

3.64e-01

</td>

<td class="gt_row gt_center">

0.37

</td>

<td class="gt_row gt_center">

196.71

</td>

<td class="gt_row gt_left">

9.51e-01

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

##### Univariate

``` r
time <- rexp(length(patients))
status <- outcome
surv.dat <- as.data.frame(cbind(time,status))
out <- uni.cox(X = gen.dat, surv.dat = surv.dat, surv.formula = Surv(time,status)~.,filter = 0.05,is.gen = T)
#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.

#> Warning: Vectorized input to `element_text()` is not officially supported.
#> Results may be unexpected or may change in future versions of ggplot2.
gt(out$tab[1:10,],rownames_to_stub = T)
```

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#hsjwkqlbsm .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
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

#hsjwkqlbsm .gt_heading {
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

#hsjwkqlbsm .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#hsjwkqlbsm .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#hsjwkqlbsm .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hsjwkqlbsm .gt_col_headings {
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

#hsjwkqlbsm .gt_col_heading {
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

#hsjwkqlbsm .gt_column_spanner_outer {
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

#hsjwkqlbsm .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#hsjwkqlbsm .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#hsjwkqlbsm .gt_column_spanner {
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

#hsjwkqlbsm .gt_group_heading {
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

#hsjwkqlbsm .gt_empty_group_heading {
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

#hsjwkqlbsm .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#hsjwkqlbsm .gt_from_md > :first-child {
  margin-top: 0;
}

#hsjwkqlbsm .gt_from_md > :last-child {
  margin-bottom: 0;
}

#hsjwkqlbsm .gt_row {
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

#hsjwkqlbsm .gt_stub {
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

#hsjwkqlbsm .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hsjwkqlbsm .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#hsjwkqlbsm .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hsjwkqlbsm .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#hsjwkqlbsm .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hsjwkqlbsm .gt_footnotes {
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

#hsjwkqlbsm .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#hsjwkqlbsm .gt_sourcenotes {
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

#hsjwkqlbsm .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#hsjwkqlbsm .gt_left {
  text-align: left;
}

#hsjwkqlbsm .gt_center {
  text-align: center;
}

#hsjwkqlbsm .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#hsjwkqlbsm .gt_font_normal {
  font-weight: normal;
}

#hsjwkqlbsm .gt_font_bold {
  font-weight: bold;
}

#hsjwkqlbsm .gt_font_italic {
  font-style: italic;
}

#hsjwkqlbsm .gt_super {
  font-size: 65%;
}

#hsjwkqlbsm .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="hsjwkqlbsm" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

Coefficient

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

Pvalue

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

MutationFrequency

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Feature

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

FDR

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left gt_stub">

MYC.Amp

</td>

<td class="gt_row gt_right">

1.4367942

</td>

<td class="gt_row gt_right">

0.003012661

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_left">

MYC.Amp

</td>

<td class="gt_row gt_right">

0.1174938

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

STK11

</td>

<td class="gt_row gt_right">

0.9302063

</td>

<td class="gt_row gt_right">

0.051591621

</td>

<td class="gt_row gt_right">

0.07

</td>

<td class="gt_row gt_left">

STK11

</td>

<td class="gt_row gt_right">

0.7208300

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

TERT

</td>

<td class="gt_row gt_right">

\-0.8621019

</td>

<td class="gt_row gt_right">

0.070297645

</td>

<td class="gt_row gt_right">

0.15

</td>

<td class="gt_row gt_left">

TERT

</td>

<td class="gt_row gt_right">

0.7208300

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

MLL2

</td>

<td class="gt_row gt_right">

0.9146774

</td>

<td class="gt_row gt_right">

0.085930560

</td>

<td class="gt_row gt_right">

0.05

</td>

<td class="gt_row gt_left">

MLL2

</td>

<td class="gt_row gt_right">

0.7208300

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

PTPRD

</td>

<td class="gt_row gt_right">

\-1.1754637

</td>

<td class="gt_row gt_right">

0.103757324

</td>

<td class="gt_row gt_right">

0.08

</td>

<td class="gt_row gt_left">

PTPRD

</td>

<td class="gt_row gt_right">

0.7208300

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

TSC1

</td>

<td class="gt_row gt_right">

\-1.6140411

</td>

<td class="gt_row gt_right">

0.110896931

</td>

<td class="gt_row gt_right">

0.05

</td>

<td class="gt_row gt_left">

TSC1

</td>

<td class="gt_row gt_right">

0.7208300

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

CDH1

</td>

<td class="gt_row gt_right">

\-1.4348192

</td>

<td class="gt_row gt_right">

0.156046843

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_left">

CDH1

</td>

<td class="gt_row gt_right">

0.7965151

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

PIK3CA

</td>

<td class="gt_row gt_right">

\-0.7835174

</td>

<td class="gt_row gt_right">

0.189345466

</td>

<td class="gt_row gt_right">

0.12

</td>

<td class="gt_row gt_left">

PIK3CA

</td>

<td class="gt_row gt_right">

0.7965151

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

EPHA5

</td>

<td class="gt_row gt_right">

\-1.2306502

</td>

<td class="gt_row gt_right">

0.223933329

</td>

<td class="gt_row gt_right">

0.05

</td>

<td class="gt_row gt_left">

EPHA5

</td>

<td class="gt_row gt_right">

0.7965151

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

FLT4

</td>

<td class="gt_row gt_right">

\-0.8701366

</td>

<td class="gt_row gt_right">

0.229947634

</td>

<td class="gt_row gt_right">

0.05

</td>

<td class="gt_row gt_left">

FLT4

</td>

<td class="gt_row gt_right">

0.7965151

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

``` r
# out$p
out$KM[[1]]
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

##### OncoCast

Axel will integrate this

#### SurvClust

Arshi will integrate this
