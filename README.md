
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
kable(gen.dat[1:10,1:10],row.names = T)
```

|                   | FLT4 | KRAS | TP53 | NF1 | ARID1A | CARD11 | ARID5B | BCOR | MLL2 | SMAD3 |
| ----------------- | ---: | ---: | ---: | --: | -----: | -----: | -----: | ---: | ---: | ----: |
| P-0010604-T01-IM5 |    0 |    0 |    1 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0002651-T01-IM3 |    0 |    0 |    1 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0000270-T01-IM3 |    0 |    0 |    1 |   0 |      1 |      1 |      1 |    0 |    0 |     0 |
| P-0002915-T01-IM3 |    0 |    0 |    0 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0011099-T01-IM5 |    0 |    0 |    0 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0000080-T01-IM3 |    0 |    0 |    1 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0001741-T01-IM3 |    0 |    1 |    1 |   0 |      1 |      0 |      0 |    0 |    0 |     0 |
| P-0003964-T01-IM3 |    0 |    1 |    0 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0003842-T01-IM5 |    0 |    0 |    0 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |
| P-0002597-T02-IM5 |    0 |    0 |    0 |   0 |      0 |      0 |      0 |    0 |    0 |     0 |

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
kable(out$fits[1:10,],row.names = T)
```

|            | Overall | 0(N=50) | 1(N=50) | OddsRatio | Pvalue   | Lower | Upper | FDR      |
| ---------- | :------ | :------ | :------ | :-------- | :------- | :---- | :---- | :------- |
| MYC.Amp    | 6%      | 12%     | 0%      | 0         | 2.67e-02 | 0     | 0.8   | 1.00e+00 |
| CDKN2B.Del | 6%      | 2%      | 10%     | 5.37      | 2.04e-01 | 0.57  | 262.3 | 1.00e+00 |
| MCL1.Amp   | 6%      | 10%     | 2%      | 0.19      | 2.04e-01 | 0     | 1.76  | 1.00e+00 |
| PIK3CA     | 12%     | 8%      | 16%     | 2.17      | 3.57e-01 | 0.53  | 10.6  | 1.00e+00 |
| FLT4       | 5%      | 2%      | 8%      | 4.21      | 3.62e-01 | 0.4   | 213.8 | 1.00e+00 |
| EPHA5      | 5%      | 8%      | 2%      | 0.24      | 3.62e-01 | 0     | 2.52  | 1.00e+00 |
| DOT1L      | 5%      | 2%      | 8%      | 4.21      | 3.62e-01 | 0.4   | 213.8 | 1.00e+00 |
| STK11      | 7%      | 10%     | 4%      | 0.38      | 4.36e-01 | 0.03  | 2.46  | 1.00e+00 |
| APC        | 7%      | 4%      | 10%     | 2.64      | 4.36e-01 | 0.41  | 29.07 | 1.00e+00 |
| MLL        | 9%      | 12%     | 6%      | 0.47      | 4.87e-01 | 0.07  | 2.37  | 1.00e+00 |

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
kable(out$tab[1:10,],row.names = T)
```

|           | Coefficient |    Pvalue | MutationFrequency | Feature   |       FDR |
| --------- | ----------: | --------: | ----------------: | :-------- | --------: |
| MLL       | \-1.4788839 | 0.0154952 |              0.09 | MLL       | 0.6043112 |
| STK11     | \-1.4567353 | 0.0519131 |              0.07 | STK11     | 0.7580415 |
| FGFR1.Amp |   1.0035559 | 0.0583109 |              0.06 | FGFR1.Amp | 0.7580415 |
| KEAP1     | \-1.0174381 | 0.1670264 |              0.05 | KEAP1     | 0.9060327 |
| NOTCH1    |   0.6087690 | 0.2071206 |              0.08 | NOTCH1    | 0.9060327 |
| DOT1L     |   0.6591677 | 0.2111284 |              0.05 | DOT1L     | 0.9060327 |
| TSC1      |   0.7504330 | 0.2131792 |              0.05 | TSC1      | 0.9060327 |
| CDH1      |   0.6047512 | 0.2520918 |              0.06 | CDH1      | 0.9060327 |
| EPHA5     | \-1.0797593 | 0.2865891 |              0.05 | EPHA5     | 0.9060327 |
| PIK3CA    |   0.4096672 | 0.2947373 |              0.12 | PIK3CA    | 0.9060327 |

``` r
# out$p
out$KM[[1]]
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

##### OncoCast

Axel will integrate this

#### SurvClust

Arshi will integrate this
