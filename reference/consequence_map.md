# Consequence Map

Data frame used as a data dictionary to recode common variant
classification types to standardized types that can be used in oncoKB
annotation.

## Usage

``` r
consequence_map
```

## Format

A data frame

- variant_classification:

  character indicating type of mutation/variant classification as it
  appears in common mutation files

- consequence_final_coding:

  final value to recode `variant_classification` column to

- consequence_final_coding_2:

  final value to recode `variant_classification` column to

- consequence_final_coding_3:

  final value to recode `variant_classification` column to

@source
<https://github.com/oncokb/oncokb-annotator/blob/a80ef0ce937c287778c36d45bf1cc8397539910c/AnnotatorCore.py#L118>
