# tbl_genomic

This function will select genes based on user inputs or on frequency
counts and then will pass the data.frame to
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html).
You can specify a `by` variable and other parameters that are accepted
by
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html).
Note the `by` variable must be merged on to onto the data before using
the `by` parameter in the function.

## Usage

``` r
tbl_genomic(
  gene_binary,
  by = NULL,
  freq_cutoff = deprecated(),
  freq_cutoff_by_gene = deprecated(),
  gene_subset = deprecated(),
  ...
)
```

## Arguments

- gene_binary:

  data.frame of genetic samples

- by:

  A variable to be passed to
  [`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)'s
  by parameter

- freq_cutoff:

  deprecated

- freq_cutoff_by_gene:

  deprecated

- gene_subset:

  deprecated

- ...:

  Additional parameters that can be passed to
  [`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html).
  To access the additional parameters you need to load `gtsummary`.

## Value

A `tbl_summary()` object

## Examples

``` r
samples <- unique(mutations$sampleId)[1:10]

gene_binary <- create_gene_binary(
  samples = samples,
  mutation = gnomeR::mutations,
  cna = gnomeR::cna,
  mut_type = "somatic_only", snp_only = FALSE,
  specify_panel = "no"
)

tbl1 <- tbl_genomic(gene_binary)

# Example wth `by` variable

gene_binary$sex <- sample( c("M", "F"), size = nrow(gene_binary), replace = TRUE)

tbl2 <- tbl_genomic(
  gene_binary = gene_binary,
  by = sex
) %>%
gtsummary::add_p() %>%
gtsummary::add_q()
#> The following errors were returned during `add_q()`:
#> ✖ For variable `AKT1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `ALK` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `APC` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `AR.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `ASXL2` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `BRCA2` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `CARD11` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `CTNNB1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `EPHB1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `FAT1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `FGFR1.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `IRS2` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `JAK1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `KDM6A.Del` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `KMT2C` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `KMT2D` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `MYC` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `MYC.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `NBN.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `NF1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `NKX3-1.Del` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `PARP1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `PDGFRA` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `PIK3R2` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `PIK3R3` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `PPP2R1A` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `PTPRS.Del` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `RECQL4.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `ROS1` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `SMAD2` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `SOX17.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `SPEN` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `SPOP` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `TP53` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `TRAF7.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
#> ✖ For variable `TSC2.Amp` (`sex`) and "estimate", "std.error", "parameter",
#>   "statistic", "conf.low", "conf.high", and "p.value" statistics: The package
#>   "broom" (>= 1.0.8) is required.
```
