# Infer variant type if not present in data

Infers variant_type from reference_allele or tumor_seq_allele data

## Usage

``` r
.infer_variant_type(mutation, names_mut_dict = names_mut_dict)
```

## Arguments

- mutation:

  mutation maf file data frame

## Value

a mutation data frame with a variant type column
