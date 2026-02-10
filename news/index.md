# Changelog

## gnomeR 1.2.2

- Cleaned up dependencies
- Removed {broom.helpers} dependency and added necessary utility
  functions as internal functions

## gnomeR 1.2.1

- More explicit separation between which are exported vs internal-only
  functions.
- Removed the sanitize_mutations(), sanitize_cna(), sanitize_fusions()
  functions and instead replaced them with smaller, more modular and
  explicit functions that do all the same tasks.
  ([\#328](https://github.com/MSKCC-Epi-Bio/gnomeR/issues/328))
- Added
  [`extract_patient_id()`](https://mskcc-epi-bio.github.io/gnomeR/reference/extract_patient_id.md)
  function to get IMPACT patient ID from sample ID
- Deprecated `freq_cutoff`, `freq_cutoff_by_gene`, and `gene_subset`
  arguments in
  [`tbl_genomic()`](https://mskcc-epi-bio.github.io/gnomeR/reference/tbl_genomic.md).
  It is now recommended that users use
  [`subset_by_frequency()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_frequency.md)
  instead before passing data to
  [`tbl_genomic()`](https://mskcc-epi-bio.github.io/gnomeR/reference/tbl_genomic.md).
- Added `other_vars` argument to
  [`subset_by_frequency()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_frequency.md),
  [`subset_by_panel()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_panel.md),
  [`summarize_by_gene()`](https://mskcc-epi-bio.github.io/gnomeR/reference/summarize_by_gene.md)
  and
  [`add_pathways()`](https://mskcc-epi-bio.github.io/gnomeR/reference/add_pathways.md)
  to allow retention of other clinical vars when using functions within
  pipeline.
- Deprecated `count_pathways_by` argument of
  [`add_pathways()`](https://mskcc-epi-bio.github.io/gnomeR/reference/add_pathways.md)
  function. Now, user must specify which specific alteration to count
  towards the pathway via the `.mut`, `.Amp`, `.Del`, `.fus` suffix
  (e.g. `custom_pathways = c('TP53.mut', 'APC.Del)`).
- Added IMPACT QA Vignette and GENIE BPC vignette
- Added `subset_by_panel` function allowing users to easily subset an
  alteration dataframe and include only genes in a specific panel
- Added IMPACT IH3 and IH4 panels to internal impact panels used for NA
  annotation in `create_gene_binary(specify_panel)`
- Added GENIE BPC alias table so users can now use
  `create_gene_binary(recode_aliases = "genie")` to check and recode
  aliases for genes in any of the GENIE BPC panels.
- Fixed bug in
  [`add_pathways()`](https://mskcc-epi-bio.github.io/gnomeR/reference/add_pathways.md)
  where `custom_pathways` wasn’t catching all types of alterations when
  `GENE.all` was used due to
  [`paste0()`](https://rdrr.io/r/base/paste.html) vectorization.
- Changed some arguments to strict matching
  ([`rlang::arg_match()`](https://rlang.r-lib.org/reference/arg_match.html))
  instead of partial matching
  ([`match.arg()`](https://rdrr.io/r/base/match.arg.html))
  (e.g. `mut_type = "s"` doesn’t work anymore and must be fully
  specified `mut_type = "somatic_only"`).
- Added unit tests for gnomeR plots/visuals
  ([\#144](https://github.com/MSKCC-Epi-Bio/gnomeR/issues/144)).
- A dictionary of old to new names for
  [`rename_columns()`](https://mskcc-epi-bio.github.io/gnomeR/reference/rename_columns.md)
  output is now an attribute of the returned object. Now messages can
  reference the original names of data columns (ex: `TumorAllele2` not
  `tumor_allele_2`) to make it more intuitive to users
  ([\#302](https://github.com/MSKCC-Epi-Bio/gnomeR/issues/302)).
- Fixed bug that wasn’t consistently filtering out germline samples
- Enhanced
  [`subset_by_frequency()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_frequency.md)
  to users to select hugo_symbols if they reach a threshold in any level
  of a variable (ex: high risk vs low risk)
  ([\#305](https://github.com/MSKCC-Epi-Bio/gnomeR/issues/305))

## gnomeR 1.2.0

- Updated color palette functionality
- Column with sample ID now returned in data frame resulting from
  [`create_gene_binary()`](https://mskcc-epi-bio.github.io/gnomeR/reference/create_gene_binary.md)
- [`summarize_by_gene()`](https://mskcc-epi-bio.github.io/gnomeR/reference/summarize_by_gene.md)
  function was changed to run faster
  ([\#259](https://github.com/MSKCC-Epi-Bio/gnomeR/issues/259))
- [`subset_by_frequency()`](https://mskcc-epi-bio.github.io/gnomeR/reference/subset_by_frequency.md)
  added to allow users to filter for specific prevalence levels of
  mutations/alterations/fusions
  ([\#270](https://github.com/MSKCC-Epi-Bio/gnomeR/issues/270))
- Removed oncoKB functionality (moving to
  [{oncokbR}](https://github.com/karissawhiting/oncokbR) package)
- Fixed bugs in
  [`tbl_genomic()`](https://mskcc-epi-bio.github.io/gnomeR/reference/tbl_genomic.md)
  and
  [`create_gene_binary()`](https://mskcc-epi-bio.github.io/gnomeR/reference/create_gene_binary.md).
- Added data processing tutorial vignette

## gnomeR 1.1.1

- Package overhauled including main binary matrix functions. This is a
  pre-v2 release made because internal workshop was taught using this
  version of functions. Final/stable versions of these functions will be
  available in v2.0.0.

## gnomeR 1.1.0

- Major changes to come. For code written pre 3/23/2022 use this release

## gnomeR 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
