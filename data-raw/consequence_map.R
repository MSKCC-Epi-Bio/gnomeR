
library(tibble)
library(dplyr)

# https://github.com/oncokb/oncokb-annotator/blob/a80ef0ce937c287778c36d45bf1cc8397539910c/AnnotatorCore.py#L118
consequence_map = tribble(
  ~variant_classification, ~consequence_final_coding, ~consequence_final_coding_2, ~consequence_final_coding_3,
  '3\'Flank',  'any', NA, NA,
  '5\'Flank',  'any', NA, NA,
  'Targeted_Region',  'inframe_deletion', 'inframe_insertion', NA,
  'COMPLEX_INDEL',  'inframe_deletion', 'inframe_insertion', NA,
  'ESSENTIAL_SPLICE_SITE',  'feature_truncation', NA, NA,
  'Exon skipping',  'inframe_deletion', NA, NA,
  'Frameshift deletion',  'frameshift_variant', NA, NA,
  'Frameshift insertion',  'frameshift_variant', NA, NA,
  'FRAMESHIFT_CODING',  'frameshift_variant', NA, NA,
  'Frame_Shift_Del',  'frameshift_variant', NA, NA,
  'Frame_Shift_Ins',  'frameshift_variant', NA, NA,
  'Fusion',  'fusion',  NA, NA,
  'Indel',  'frameshift_variant', 'inframe_deletion', 'inframe_insertion',
  'In_Frame_Del',  'inframe_deletion', NA, NA,
  'In_Frame_Ins',  'inframe_insertion', NA, NA,
  'Missense',  'missense_variant', NA, NA,
  'Missense_Mutation',  'missense_variant', NA, NA,
  'Nonsense_Mutation',  'stop_gained', NA, NA,
  'Nonstop_Mutation',  'stop_lost', NA, NA,
  'Splice_Site',  'splice_region_variant', NA, NA,
  'Splice_Site_Del',  'splice_region_variant', NA, NA,
  'Splice_Site_SNP',  'splice_region_variant', NA, NA,
  'splicing',  'splice_region_variant', NA, NA,
  'Translation_Start_Site',  'start_lost', NA, NA,
  'vIII deletion',  'any', NA, NA,

  # Karissa Added,
  'Splice_Region',  'splice_region_variant', NA, NA,
  'Intron', 'intron_variant', NA, NA
)

consequence_map <- consequence_map %>%
  mutate(across(everything(), ~tolower(.x)))

usethis::use_data(consequence_map, overwrite = TRUE)
