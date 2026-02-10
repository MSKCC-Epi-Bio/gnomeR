# Sum Alterations in a Pathway

Sum Alterations in a Pathway

## Usage

``` r
.sum_alts_in_pathway(
  gene_binary,
  pathway_list_item,
  pathway_name,
  count_pathways_by
)
```

## Arguments

- gene_binary:

  a binary matrix (see `gene_binary()`)

- pathway_list_item:

  a named list of length 1 with pathway name as name and vector of genes
  as first and only item in list

- pathway_name:

  name of pathway

## Value

a dataframe of 1 column of 0/1s indicating pathway activated yes/no

## Examples

``` r
gene_binary <- create_gene_binary(mutation = gnomeR::mutations, cna = gnomeR::cna,
fusion = gnomeR::sv)
#> ! `samples` argument is `NULL`. We will infer your cohort inclusion and resulting data frame will include all samples with at least one alteration in mutation, fusion or cna data frames
#> ! 7 mutations have `NA` or blank in the mutationStatus column instead of 'SOMATIC' or 'GERMLINE'. These were assumed to be 'SOMATIC' and were retained in the resulting binary matrix.
x <- .sum_alts_in_pathway(gene_binary,
 pathway_list_item = gnomeR::pathways[1],
  pathway_name = names(gnomeR::pathways[1]),
    count_pathways_by = "alteration")
```
