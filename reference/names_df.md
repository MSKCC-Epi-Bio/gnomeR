# Data Frame of Column Names

Data frame of accepted data names for standard genomic files. This
serves as a dictionary to help disambiguate raw column names from user
entered mutation, CNA or structural variant data

## Usage

``` r
names_df
```

## Format

A data frame

- maf_column_name:

  data field names as they appear in common MAF file

- api_column_name:

  data field names as they appear in common cBioPortal API retrieved
  files

- mutation_input:

  does this field appear in mutation files?

- fusion_input:

  does this field appear in mutation/sv files?

- cna_input:

  does this field appear in CNA files?

- definition:

  variable definition

- notes:

  data notes

- sc_maf_column_name:

  snake case version of `maf_column_name`

- sc_api_column_name:

  snake case version of `api_column_name`

- internal_column_name:

  name used for each field for all internal processing functions
