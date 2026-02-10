# An example IMPACT cBioPortal mutation data set in API format

This set was created from a random sample of 200 patients from publicly
available prostate cancer data from cBioPortal. The file is in API
format.

## Usage

``` r
sv
```

## Format

A data frame with structural variants from Abida et al. JCO Precis Oncol
2017. Retrieved from cBioPortal.There are 94 observations and 29
variables.

A data frame with 94 rows and 44 variables:

- `uniqueSampleKey`:

  character COLUMN_DESCRIPTION

- `uniquePatientKey`:

  character COLUMN_DESCRIPTION

- `molecularProfileId`:

  character COLUMN_DESCRIPTION

- `sampleId`:

  MSKCC Sample ID

- `patientId`:

  Patient ID

- `studyId`:

  Indicator for Abida et al. 2017 study

- `site1EntrezGeneId`:

  integer COLUMN_DESCRIPTION

- `site1HugoSymbol`:

  Character w/ 31 levels, Column containing HUGO symbols genes for first
  site of fusion

- `site1EnsemblTranscriptId`:

  character COLUMN_DESCRIPTION

- `site1Chromosome`:

  character COLUMN_DESCRIPTION

- `site1Position`:

  integer COLUMN_DESCRIPTION

- `site1Contig`:

  character COLUMN_DESCRIPTION

- `site1Region`:

  character COLUMN_DESCRIPTION

- `site1RegionNumber`:

  integer COLUMN_DESCRIPTION

- `site1Description`:

  character COLUMN_DESCRIPTION

- `site2EntrezGeneId`:

  integer COLUMN_DESCRIPTION

- `site2HugoSymbol`:

  Character w/ 21 levels, Column containing all HUGO symbols genes for
  second site of fusion

- `site2EnsemblTranscriptId`:

  character COLUMN_DESCRIPTION

- `site2Chromosome`:

  character COLUMN_DESCRIPTION

- `site2Position`:

  integer COLUMN_DESCRIPTION

- `site2Contig`:

  character COLUMN_DESCRIPTION

- `site2Region`:

  character COLUMN_DESCRIPTION

- `site2RegionNumber`:

  integer COLUMN_DESCRIPTION

- `site2Description`:

  character COLUMN_DESCRIPTION

- `site2EffectOnFrame`:

  character COLUMN_DESCRIPTION

- `ncbiBuild`:

  character COLUMN_DESCRIPTION

- `dnaSupport`:

  Factor, all are `yes` in this data

- `rnaSupport`:

  Factor, all are `unknown` in this data

- `normalReadCount`:

  integer COLUMN_DESCRIPTION

- `tumorReadCount`:

  integer COLUMN_DESCRIPTION

- `normalVariantCount`:

  integer COLUMN_DESCRIPTION

- `tumorVariantCount`:

  integer COLUMN_DESCRIPTION

- `normalPairedEndReadCount`:

  integer COLUMN_DESCRIPTION

- `tumorPairedEndReadCount`:

  integer COLUMN_DESCRIPTION

- `normalSplitReadCount`:

  integer COLUMN_DESCRIPTION

- `tumorSplitReadCount`:

  integer COLUMN_DESCRIPTION

- `annotation`:

  character COLUMN_DESCRIPTION

- `breakpointType`:

  character COLUMN_DESCRIPTION

- `connectionType`:

  character COLUMN_DESCRIPTION

- `eventInfo`:

  character COLUMN_DESCRIPTION

- `variantClass`:

  character COLUMN_DESCRIPTION

- `length`:

  integer COLUMN_DESCRIPTION

- `comments`:

  character COLUMN_DESCRIPTION

- `svStatus`:

  character COLUMN_DESCRIPTION

## Source

<https://www.cbioportal.org/study/summary?id=prad_mskcc_2017>
