# Simplify binary matrix to one column per patient that counts any alteration type across all samples as 1

This will reduce the number of columns in your binary matrix, and the
resulting data frame will have only 1 col per gene, as opposed to
separate columns for mutation/cna/fusion.

## Usage

``` r
summarize_by_patient(gene_binary, other_vars = NULL)
```

## Arguments

- gene_binary:

  a 0/1 matrix of gene alterations

- other_vars:

  One or more column names (quoted or unquoted) in data to be retained
  in resulting data frame. Default is NULL.

## Value

a binary matrix with a row for each sample and one column per gene

## Details

Note that if samples to the same patient were sequenced on different
panels, any indication of an alteration is counted as an alteration, but
the absence of an alteration is only defined when all sequencing panels
included the gene and indicated that it was not altered.

## Examples

``` r
samples <- unique(gnomeR::mutations$sampleId)[1:10]
gene_binary <- create_gene_binary(
  samples = samples, mutation = mutations, cna = cna,
  mut_type = "somatic_only",
  include_silent = FALSE,
  specify_panel = "IMPACT341")

gene_binary$patient_id = extract_patient_id(gene_binary$sample_id)

summarize_by_patient(gene_binary)
#> # A tibble: 9 × 36
#>   patient_id PARP1  AKT1   ALK   APC BRCA2 CTNNB1 EPHB1  FAT1  JAK1 SMAD2   NF1
#>   <chr>      <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 P-0001128      1     1     0     0     1      1     0     0     0     0     0
#> 2 P-0001859      1     0     0     1     0      0     1     0     0     0     0
#> 3 P-0001895      1     0     1     0     0      0     0     0     0     0     0
#> 4 P-0001845      0     1     0     0     0      0     0     0     0     0     1
#> 5 P-0005570      0     1     0     0     0      0     0     0     0     0     0
#> 6 P-0001768      0     0     1     0     0      0     0     0     0     0     0
#> 7 P-0004508      0     0     1     0     0      0     0     0     0     0     0
#> 8 P-0002984      0     0     1     0     0      0     0     0     1     0     0
#> 9 P-0000964      0     0     0     1     1      0     0     1     0     1     0
#> # ℹ 24 more variables: PDGFRA <dbl>, PIK3R2 <dbl>, PPP2R1A <dbl>, ROS1 <dbl>,
#> #   TP53 <dbl>, KMT2D <dbl>, SPOP <dbl>, PIK3R3 <dbl>, IRS2 <dbl>, SPEN <dbl>,
#> #   ASXL2 <dbl>, KMT2C <dbl>, CARD11 <dbl>, PTPRS <dbl>, KDM6A <dbl>,
#> #   `NKX3-1` <dbl>, AR <dbl>, TRAF7 <dbl>, TSC2 <dbl>, FGFR1 <dbl>,
#> #   SOX17 <dbl>, RECQL4 <dbl>, NBN <dbl>, MYC <dbl>
```
