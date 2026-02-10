# An example IMPACT cBioPortal mutation data set in API format

This set contains a random sample of 200 patients from publicly
available prostate cancer data from cBioPortal. The file is in API
format.

## Usage

``` r
mutations
```

## Format

A data frame with mutations from Abida et al. JCO Precis Oncol 2017.
Retrieved from cBioPortal.There are 725 observations and 29 variables.

- hugoGeneSymbol:

  Character w/ 324 levels, Column containing all HUGO symbols genes

- entrezGeneId:

  Entrez Gene ID

- sampleId:

  MSKCC Sample ID

- patientId:

  Patient ID

- studyId:

  Indicator for Abida et al. 2017 study

- center:

  Cancer Center ID

- mutationStatus:

  Somatic or germ-line mutation status

- variantType:

  Mutation variant type

- chr:

  Chromosome mutation observed on

- endPosition:

  End Position

- fisValue:

  fisValue

- functionalImpactScore:

  functionalImpactScore

- keyword:

  keyword

- linkMsa:

  linkMsa

- linkPdb:

  linkPdb

- linkXvar:

  linkXvar

- molecularProfileId:

  molecularProfileId

- mutationType:

  mutationType

- ncbiBuild:

  ncbiBuild

- proteinChange:

  proteinChange

- proteinPosEnd:

  proteinPosEnd

- proteinPosStart:

  proteinPosStart

- referenceAllele:

  referenceAllele

- refseqMrnaId:

  refseqMrnaId

- startPosition:

  startPosition

- uniquePatientKey:

  uniquePatientKey

- uniqueSampleKey:

  uniqueSampleKey

- validationStatus:

  validationStatus

- variantAllele:

  variantAllele

## Source

<https://www.cbioportal.org/study/summary?id=prad_mskcc_2017>
