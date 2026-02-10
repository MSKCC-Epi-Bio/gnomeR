# Annotate Missing Gene Values According to Specific Panels

Annotate Missing Gene Values According to Specific Panels

## Usage

``` r
annotate_any_panel(sample_panel_pair, gene_binary)
```

## Arguments

- sample_panel_pair:

  a data frame of `sample_id`-`panel_id` pairs specifying panels to use
  for annotation of each sample

- gene_binary:

  a binary matrix of 0/1 indicating alteration yes/no for each sample

## Value

a gene_binary annotated for missingness
