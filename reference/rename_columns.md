# Rename columns from API results to work with gnomeR functions

Will return a named vector of internal column names as values and
original data set names as names as an attribute
(`attr(x, "names_dict")`)

## Usage

``` r
rename_columns(df_to_check)
```

## Arguments

- df_to_check:

  A data frame to check and recode names as needed

## Value

a renamed data frame

## Examples

``` r
rename_columns(df_to_check = gnomeR::mutations)
#> # A tibble: 725 × 29
#>    hugo_symbol entrez_gene_id uniqueSampleKey                   uniquePatientKey
#>    <chr>                <int> <chr>                             <chr>           
#>  1 PARP1                  142 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9t… UC0wMDAxMTI4OnB…
#>  2 PARP1                  142 UC0wMDAxODU5LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODU5OnB…
#>  3 PARP1                  142 UC0wMDAxODk1LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODk1OnB…
#>  4 AKT1                   207 UC0wMDAxMTI4LVQwMS1JTTM6cHJhZF9t… UC0wMDAxMTI4OnB…
#>  5 AKT1                   207 UC0wMDAxODQ1LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODQ1OnB…
#>  6 AKT1                   207 UC0wMDA1NTcwLVQwMS1JTTU6cHJhZF9t… UC0wMDA1NTcwOnB…
#>  7 ALK                    238 UC0wMDAxNzY4LVQwMS1JTTM6cHJhZF9t… UC0wMDAxNzY4OnB…
#>  8 ALK                    238 UC0wMDA0NTA4LVQwMS1JTTU6cHJhZF9t… UC0wMDA0NTA4OnB…
#>  9 ALK                    238 UC0wMDAxODk1LVQwMS1JTTM6cHJhZF9t… UC0wMDAxODk1OnB…
#> 10 ALK                    238 UC0wMDAyOTg0LVQwMS1JTTM6cHJhZF9t… UC0wMDAyOTg0OnB…
#> # ℹ 715 more rows
#> # ℹ 25 more variables: molecular_profile_id <chr>, sample_id <chr>,
#> #   patient_id <chr>, study_id <chr>, center <chr>, mutation_status <chr>,
#> #   validation_status <chr>, start_position <int>, end_position <int>,
#> #   reference_allele <chr>, hgv_sp_short <chr>, variant_classification <chr>,
#> #   functionalImpactScore <chr>, fisValue <dbl>, linkXvar <chr>, linkPdb <chr>,
#> #   linkMsa <chr>, ncbi_build <chr>, variant_type <chr>, keyword <chr>, …
x <- rename_columns(df_to_check = gnomeR::sv)
attr(x, "names_dict")
#>                    sample_id                   patient_id 
#>                   "sampleId"                  "patientId" 
#>                     study_id         molecular_profile_id 
#>                    "studyId"         "molecularProfileId" 
#>                   ncbi_build           site_1_hugo_symbol 
#>                  "ncbiBuild"            "site1HugoSymbol" 
#>           site_2_hugo_symbol       site_2_effect_on_frame 
#>            "site2HugoSymbol"         "site2EffectOnFrame" 
#>                    sv_status                  dna_support 
#>                   "svStatus"                 "dnaSupport" 
#>                  rna_support                     comments 
#>                 "rnaSupport"                   "comments" 
#>                   event_info        site_1_entrez_gene_id 
#>                  "eventInfo"          "site1EntrezGeneId" 
#> site_1_ensembl_transcript_id            site_1_chromosome 
#>   "site1EnsemblTranscriptId"            "site1Chromosome" 
#>              site_1_position                site_1_contig 
#>              "site1Position"                "site1Contig" 
#>                site_1_region         site_1_region_number 
#>                "site1Region"          "site1RegionNumber" 
#>           site_1_description        site_2_entrez_gene_id 
#>           "site1Description"          "site2EntrezGeneId" 
#> site_2_ensembl_transcript_id            site_2_chromosome 
#>   "site2EnsemblTranscriptId"            "site2Chromosome" 
#>              site_2_position                site_2_contig 
#>              "site2Position"                "site2Contig" 
#>                site_2_region         site_2_region_number 
#>                "site2Region"          "site2RegionNumber" 
#>           site_2_description                        class 
#>           "site2Description"               "variantClass" 
```
