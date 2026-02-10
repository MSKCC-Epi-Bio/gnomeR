# Pivot CNA from maf (long) version to wide version

Pivot CNA from maf (long) version to wide version

## Usage

``` r
pivot_cna_wider(cna)
```

## Arguments

- cna:

  a cna dataframe in maf (long) format

## Value

a dataframe of reformatted CNA alteration (in wide format)

## Examples

``` r
cna_long <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c("AMPLIFICATION","AMPLIFICATION",
                   "AMPLIFICATION","AMPLIFICATION","DELETION"))

cna <- pivot_cna_wider(cna_long)

 cna_long <- data.frame(
    sampleId = c("P-0001276-T01-IM3","P-0001276-T01-IM3",
                 "P-0005436-T01-IM3",
                 "P-0001276-T01-IM3","P-0003333-T01-IM3"),
    Hugo_Symbol = c("MLL2","KMT2D","HIST1H2BD",
                    "HIST1H3B","KDR"),
    alteration = c(2, 2, -1, 1, -2))

cna <- pivot_cna_wider(cna_long)
```
