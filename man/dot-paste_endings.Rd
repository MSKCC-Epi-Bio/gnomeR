% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{.paste_endings}
\alias{.paste_endings}
\title{Add descriptive endings to hugo symbol names that do not have one already}
\usage{
.paste_endings(names, ending = NULL)
}
\arguments{
\item{names}{hugo symbols to check}

\item{ending}{character ending to add to hugo symbol names without descriptive endings.
The default is ".mut". If interested in any type of alteration, use ".any".}
}
\value{
a vector of hugo symbols where each entry has a descriptive ending
from the following list: ".Amp", ".Del", ".fus", ".cna", ".mut".
}
\description{
Add descriptive endings to hugo symbol names that do not have one already
}
\keyword{internal}
