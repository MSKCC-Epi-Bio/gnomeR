# gnomeR (development version)

- Deprecated `freq_cutoff`, `freq_cutoff_by_gene`, and `gene_subset` arguments in `tbl_genomic()`. It is now recommended that users use `subset_by_frequency()` instead before passing data to `tbl_genomic()`.
- Added `other_vars` argument to `subset_by_frequency()` to allow retention of other clinical vars when using function within pipeline.
- Deprecated `count_pathways_by` argument of `add_pathways()` function. Now, user must specify which specific alteration to count towards the pathway via the `.mut`, `.Amp`, `.Del`, `.fus` suffix (e.g. `custom_pathways = c('TP53.mut', 'APC.Del)`). 
- Added IMPACT QA Vignette and GENIE BPC vignette

# gnomeR 1.2.0

* Updated color palette functionality
* Column with sample ID now returned in data frame resulting from `create_gene_binary()`
* `summarize_by_gene()` function was changed to run faster (#259)
* `subset_by_frequency()` added to allow users to filter for specific prevalence levels of mutations/alterations/fusions (#270)
* Removed oncoKB functionality (moving to [{oncokbR}](https://github.com/karissawhiting/oncokbR) package)
* Fixed bugs in `tbl_genomic()` and `create_gene_binary()`. 
* Added data processing tutorial vignette 

# gnomeR 1.1.1

* Package overhauled including main binary matrix functions. This is a pre-v2 release made because internal workshop was taught using this version of functions. Final/stable versions of these functions will be available in v2.0.0.

# gnomeR 1.1.0

* Major changes to come. For code written pre 3/23/2022 use this release


# gnomeR 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
