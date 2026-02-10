# Subset a Binary Matrix By Genes Available on Specified Panel

Subset a Binary Matrix By Genes Available on Specified Panel

## Usage

``` r
subset_by_panel(gene_binary, panel_id = NULL, other_vars = NULL)
```

## Arguments

- gene_binary:

  A data frame with a row for each sample and column for each
  alteration. Data frame must have a `sample_id` column and columns for
  each alteration with values of 0, 1 or NA.

- panel_id:

  A character string or vector of the specified panel to subset the
  genes (see
  [`gnomeR::gene_panels`](https://mskcc-epi-bio.github.io/gnomeR/reference/gene_panels.md)
  for available panels)

- other_vars:

  One or more column names (quoted or unquoted) in data to be retained
  in resulting data frame. Default is NULL.

## Value

a data frame with a `sample_id` column and columns for alterations on
genes that were sequenced on the specified panel.

## Author

Jessica Lavery

## Examples

``` r
samples <- unique(gnomeR::mutations$sampleId)
gene_binary <- create_gene_binary(
  samples = samples, mutation = mutations, cna = cna,
  mut_type = "somatic_only",
  include_silent = FALSE,
  specify_panel = "impact"
)
subset_by_panel(gene_binary = gene_binary, panel_id = "IMPACT468")
#> ℹ There are 244 genes on the IMPACT468 panel that are not altered for any patients; columns for those genes are not included in the resulting data frame. To see the names of the genes that are not mutated for any patients in the sample, To view and compare genes in each panel, see `unnest(gnomeR::gene_panels, cols = c('genes_in_panel'))`
#> # A tibble: 200 × 228
#>    sample_id      AKT1  AKT3   ALK ANKRD11   APC    AR  ARAF ARID1A ARID1B ARID2
#>    <chr>         <dbl> <dbl> <dbl>   <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>
#>  1 P-0001128-T0…     1     0     0      NA     0     0     0      0      0     0
#>  2 P-0001859-T0…     0     0     0      NA     1     0     0      0      0     0
#>  3 P-0001895-T0…     0     0     1      NA     0     0     0      0      0     0
#>  4 P-0001845-T0…     1     0     0      NA     0     0     0      0      0     0
#>  5 P-0001768-T0…     0     0     1      NA     0     0     0      0      0     0
#>  6 P-0002984-T0…     0     0     1      NA     0     0     0      0      0     0
#>  7 P-0000964-T0…     0     0     0      NA     1     0     0      0      0     0
#>  8 P-0000964-T0…     0     0     0      NA     1     0     0      0      0     0
#>  9 P-0000610-T0…     0     0     0      NA     1     0     0      0      0     0
#> 10 P-0001247-T0…     0     0     0      NA     1     0     0      0      0     0
#> # ℹ 190 more rows
#> # ℹ 217 more variables: ARID5B <dbl>, ASXL1 <dbl>, ASXL2 <dbl>, ATM <dbl>,
#> #   ATR <dbl>, ATRX <dbl>, AURKB <dbl>, AXIN1 <dbl>, AXIN2 <dbl>, BAP1 <dbl>,
#> #   BARD1 <dbl>, BCL6 <dbl>, BCOR <dbl>, BLM <dbl>, BMPR1A <dbl>, BRAF <dbl>,
#> #   BRCA1 <dbl>, BRCA2 <dbl>, BRIP1 <dbl>, CARD11 <dbl>, CBL <dbl>,
#> #   CCND1 <dbl>, CD276 <dbl>, CD79B <dbl>, CDC73 <dbl>, CDH1 <dbl>,
#> #   CDK12 <dbl>, CDK4 <dbl>, CDK8 <dbl>, CDKN1A <dbl>, CDKN2C <dbl>, …

p_genes <- tidyr::unnest(gnomeR::gene_panels, cols = c("genes_in_panel"))
p_genes <- p_genes[p_genes$gene_panel == 'IMPACT300', ]
setdiff(p_genes$genes_in_panel, names(gene_binary))
#>   [1] "ABL1"     "ABL2"     "AKT2"     "ALOX12B"  "AMER1"    "ARHGAP26"
#>   [7] "AURKA"    "BCL2L1"   "BCL2L11"  "BIRC2"    "BUB1B"    "CBLB"    
#>  [13] "CBLC"     "CCND2"    "CCND3"    "CCNE1"    "CDC42EP2" "CDH11"   
#>  [19] "CDK6"     "CDKN2A"   "CDKN2B"   "CEBPA"    "CHEK1"    "CRKL"    
#>  [25] "CRLF2"    "CYLD"     "DIS3"     "DNMT3A"   "DNMT3B"   "E2F3"    
#>  [31] "EIF4EBP1" "EPHA10"   "EPHA2"    "EPHA4"    "EPHA6"    "EPHA7"   
#>  [37] "EPHA8"    "EPHB2"    "EPHB3"    "EPHB4"    "EPHB6"    "EZH2"    
#>  [43] "FAS"      "FAT4"     "FBXO11"   "FBXW7"    "FGFR3"    "FKBP1A"  
#>  [49] "FLT4"     "GATA1"    "GATA2"    "GLI3"     "GNA11"    "GNAS"    
#>  [55] "GOLPH3"   "GRM3"     "GSK3B"    "HDAC2"    "HIF1A"    "HLA-A"   
#>  [61] "HMGA2"    "HSP90AA1" "IGFBP7"   "IKBKE"    "IL7R"     "INSR"    
#>  [67] "IRF4"     "JUN"      "KCNJ5"    "KDR"      "KEAP1"    "KIT"     
#>  [73] "KLF6"     "LDHA"     "LGR6"     "MAGI2"    "MAP2K4"   "MAP3K1"  
#>  [79] "MAP3K8"   "MCL1"     "MLST8"    "MPL"      "MYB"      "MYCL"    
#>  [85] "MYCN"     "MYD88"    "NCOA2"    "NF2"      "NFE2L2"   "NFKB1"   
#>  [91] "NFKB2"    "NOTCH2"   "NRAS"     "PAK5"     "PBRM1"    "PHOX2B"  
#>  [97] "PKM"      "PNRC1"    "POLE"     "PPP6C"    "PREX2"    "PRKAA2"  
#> [103] "PRKAR1A"  "PRKCI"    "PRKN"     "PTPN11"   "RAC1"     "RAD50"   
#> [109] "RAF1"     "REL"      "RICTOR"   "ROR2"     "RPS6KB1"  "SMARCA4" 
#> [115] "SMARCB1"  "SOCS1"    "SRC"      "SRSF2"    "SUFU"     "SYK"     
#> [121] "TBK1"     "TEK"      "TERT"     "TNFRSF14" "TSC1"     "TSHR"    
#> [127] "WAS"      "WNK1"     "YES1"     "ZRSR2"   
```
