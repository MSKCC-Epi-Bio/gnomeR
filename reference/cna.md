# An example IMPACT cBioPortal mutation data set in API format

This set was created from a random sample of 200 patients from publicly
available prostate cancer data from cBioPortal. The file is in API
format.

## Usage

``` r
cna
```

## Format

A data frame with copy number alterations (CNA) from Abida et al. JCO
Precis Oncol 2017.Retrieved from cBioPortal.There are 475 observations
and 29 variables.

- hugoGeneSymbol:

  Character w/ 324 levels, Column containing all HUGO symbols genes

- entrezGeneId:

  Entrez Gene ID

- molecularProfileId:

  Molecular Profile ID for data set

- sampleId:

  MSKCC Sample ID

- patientId:

  Patient ID

- studyId:

  Indicator for Abida et al. 2017 study

- alteration:

  Factor, Type of CNA

- uniqueSampleKey:

  character COLUMN_DESCRIPTION

- uniquePatientKey:

  character COLUMN_DESCRIPTION

## Source

<https://www.cbioportal.org/study/summary?id=prad_mskcc_2017>
