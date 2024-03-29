% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sv}
\alias{sv}
\title{An example IMPACT cBioPortal mutation data set in API format}
\format{
A data frame with structural variants from Abida et al. JCO Precis Oncol 2017.
Retrieved from cBioPortal.There are 94 observations and 29 variables.

A data frame with 94 rows and 44 variables:
\describe{
\item{\code{uniqueSampleKey}}{character COLUMN_DESCRIPTION}
\item{\code{uniquePatientKey}}{character COLUMN_DESCRIPTION}
\item{\code{molecularProfileId}}{character COLUMN_DESCRIPTION}
\item{\code{sampleId}}{MSKCC Sample ID}
\item{\code{patientId}}{Patient ID}
\item{\code{studyId}}{Indicator for Abida et al. 2017 study}
\item{\code{site1EntrezGeneId}}{integer COLUMN_DESCRIPTION}
\item{\code{site1HugoSymbol}}{Character w/ 31 levels,
Column containing HUGO symbols genes for first site of fusion}
\item{\code{site1EnsemblTranscriptId}}{character COLUMN_DESCRIPTION}
\item{\code{site1Chromosome}}{character COLUMN_DESCRIPTION}
\item{\code{site1Position}}{integer COLUMN_DESCRIPTION}
\item{\code{site1Contig}}{character COLUMN_DESCRIPTION}
\item{\code{site1Region}}{character COLUMN_DESCRIPTION}
\item{\code{site1RegionNumber}}{integer COLUMN_DESCRIPTION}
\item{\code{site1Description}}{character COLUMN_DESCRIPTION}
\item{\code{site2EntrezGeneId}}{integer COLUMN_DESCRIPTION}
\item{\code{site2HugoSymbol}}{Character w/ 21 levels,
Column containing all HUGO symbols genes for second site of fusion}
\item{\code{site2EnsemblTranscriptId}}{character COLUMN_DESCRIPTION}
\item{\code{site2Chromosome}}{character COLUMN_DESCRIPTION}
\item{\code{site2Position}}{integer COLUMN_DESCRIPTION}
\item{\code{site2Contig}}{character COLUMN_DESCRIPTION}
\item{\code{site2Region}}{character COLUMN_DESCRIPTION}
\item{\code{site2RegionNumber}}{integer COLUMN_DESCRIPTION}
\item{\code{site2Description}}{character COLUMN_DESCRIPTION}
\item{\code{site2EffectOnFrame}}{character COLUMN_DESCRIPTION}
\item{\code{ncbiBuild}}{character COLUMN_DESCRIPTION}
\item{\code{dnaSupport}}{Factor, all are \code{yes} in this data}
\item{\code{rnaSupport}}{Factor, all are \code{unknown} in this data}
\item{\code{normalReadCount}}{integer COLUMN_DESCRIPTION}
\item{\code{tumorReadCount}}{integer COLUMN_DESCRIPTION}
\item{\code{normalVariantCount}}{integer COLUMN_DESCRIPTION}
\item{\code{tumorVariantCount}}{integer COLUMN_DESCRIPTION}
\item{\code{normalPairedEndReadCount}}{integer COLUMN_DESCRIPTION}
\item{\code{tumorPairedEndReadCount}}{integer COLUMN_DESCRIPTION}
\item{\code{normalSplitReadCount}}{integer COLUMN_DESCRIPTION}
\item{\code{tumorSplitReadCount}}{integer COLUMN_DESCRIPTION}
\item{\code{annotation}}{character COLUMN_DESCRIPTION}
\item{\code{breakpointType}}{character COLUMN_DESCRIPTION}
\item{\code{connectionType}}{character COLUMN_DESCRIPTION}
\item{\code{eventInfo}}{character COLUMN_DESCRIPTION}
\item{\code{variantClass}}{character COLUMN_DESCRIPTION}
\item{\code{length}}{integer COLUMN_DESCRIPTION}
\item{\code{comments}}{character COLUMN_DESCRIPTION}
\item{\code{svStatus}}{character COLUMN_DESCRIPTION}
}
}
\source{
\url{https://www.cbioportal.org/study/summary?id=prad_mskcc_2017}
}
\usage{
sv
}
\description{
This set was created from a random sample of 200 patients from
publicly available prostate cancer data from cBioPortal. The file
is in API format.
}
\keyword{datasets}
