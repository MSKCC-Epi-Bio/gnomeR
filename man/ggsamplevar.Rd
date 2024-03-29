% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting-functions.R
\name{ggsamplevar}
\alias{ggsamplevar}
\title{#' Utility Function to Extract SNV
#'
#' @param x string
#' @param n number of characters from right
#'
#' @return string
#' @noRd
#' @examples
#' substrRight("Hello", 2)
#'
substrRight <- function(x, n) {
x <- as.character(x)
substr(x, nchar(x) - n + 1, nchar(x))
}
Histogram of Variants Per Sample Colored By Variant Classification}
\usage{
ggsamplevar(mutation)
}
\arguments{
\item{mutation}{Raw mutation dataframe containing alteration data}
}
\value{
Histogram of counts of variants per tumor sample
}
\description{
#' Utility Function to Extract SNV
#'
#' @param x string
#' @param n number of characters from right
#'
#' @return string
#' @noRd
#' @examples
#' substrRight("Hello", 2)
#'
substrRight <- function(x, n) {
x <- as.character(x)
substr(x, nchar(x) - n + 1, nchar(x))
}
Histogram of Variants Per Sample Colored By Variant Classification
}
\examples{
ggsamplevar(gnomeR::mutations)

}
