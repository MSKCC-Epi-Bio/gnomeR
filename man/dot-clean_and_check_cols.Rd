% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-gene-binary.R
\name{.clean_and_check_cols}
\alias{.clean_and_check_cols}
\title{Checks genomic input file columns to ensure column names are correct}
\usage{
.clean_and_check_cols(
  df_to_check,
  required_cols = c("sample_id", "hugo_symbol")
)
}
\arguments{
\item{df_to_check}{Raw maf dataframe containing alteration data}

\item{required_cols}{A character specifying names of columns to check}
}
\value{
a corrected maf file or an error if problems with maf
}
\description{
Checks genomic input file columns to ensure column names are correct
}
\examples{
gnomeR:::.clean_and_check_cols(df_to_check = gnomeR::mutations)

}
\keyword{internal}
