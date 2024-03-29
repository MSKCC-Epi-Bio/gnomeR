% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create-gene-binary.R
\name{.cna_gene_binary}
\alias{.cna_gene_binary}
\title{Make Binary Matrix From CNA data frame}
\usage{
.cna_gene_binary(cna, samples, specify_panel, high_level_cna_only)
}
\arguments{
\item{cna}{A data frame of copy number alterations. If inputed the outcome will be added to the matrix with columns ending in ".del" and ".amp".
Default is NULL.}

\item{samples}{a character vector specifying which samples should be included in the resulting data frame.
Default is NULL is which case all samples with an alteration in the mutation, cna or fusions file will be used. If you specify
a vector of samples that contain samples not in any of the passed genomic data frames, 0's (or NAs when appropriate if specifying a panel) will be
returned for every column of that patient row.}

\item{specify_panel}{Default is \code{"no"} where no panel annotation is done. Otherwise pass a character vector of length 1 with a
panel id (see \code{gnomeR::gene_panels} for available panels), or \code{"impact"} for automated IMPACT annotation.
Alternatively, you may pass a data frame of \code{sample_id}-\code{panel_id} pairs specifying panels for each sample for
which to insert NAs indicating genes not tested. See below for details.}

\item{high_level_cna_only}{If TRUE, only deep deletions (-2, -1.5) or high level amplifications (2) will be counted as events
in the binary matrix. Gains (1) and losses (1) will be ignored. Default is \code{FALSE} where all CNA events are counted.}
}
\value{
a data frame
}
\description{
Make Binary Matrix From CNA data frame
}
\keyword{internal}
