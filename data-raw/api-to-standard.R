library(cbioportalR)
library(dplyr)

rename_df <- tibble::tribble(
                           ~maf_col,      ~api_col,
                     "Hugo_Symbol",     "hugoGeneSymbol",
                  "Entrez_Gene_Id",       "entrezGeneId",
                          "Center",             "center",
                      "NCBI_Build",          "ncbiBuild",
                      "Chromosome",                "chr",
                  "Start_Position",      "startPosition",
                    "End_Position",        "endPosition",
          "Variant_Classification",       "mutationType",
                   "Variant_Type",        "variantType",
               "Reference_Allele",    "referenceAllele",
           "Tumor_Sample_Barcode",           "sampleId",
              "Validation_Status",   "validationStatus",
                "Mutation_Status",     "mutationStatus",
                    "HGVSp_Short",      "proteinChange",
                         "Allele",      "variantAllele"
  ) %>%
  mutate(api_snakecase = janitor::make_clean_names(api_col))



usethis::use_data(rename_df, overwrite = TRUE)
