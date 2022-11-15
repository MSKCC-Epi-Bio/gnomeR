#' IMPACT Gene Meta Data
#'
#' Dataframe labeling all genes found in IMPACT along with their corresponding
#' platform and Entrez ID.
#'
#' @format A data frame with 469 observations and 5 variables
#' \describe{
#'   \item{hugo_symbol}{Column containing all HUGO symbols genes included in IMPACT}
#'   \item{platform_341}{Character, indicates whether each gene was included in
#'   IMPACT platform 341. Options are `included` and `not included`}
#'   \item{platform_410}{Character, indicates whether each gene was included in
#'   IMPACT platform 410. Options are `included` and `not included`}
#'   \item{platform_468}{Character, indicates whether each gene was included in
#'   IMPACT platform 468. Options are `included` and `not included`}
#'   \item{alias}{List of EntrezIDs corresponding to the value in the hugo_symbol column and alias list.}
#'   \item{entrez_id}{List of known aliases corresponding to the value in the hugo_symbol column.}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"impact_gene_info"

#' Public Gene Panels on cBioPortal
#'
#' Data frame of all gene panels available in public cBioPortal
#'
#' @format A nested data frame
#' \describe{
#'   \item{gene_panels}{Gene panel ID}
#'   \item{genes_in_panel}{List column of Hugo symbols of all genes in each panel}
#'   \item{entrez_ids_in_panel}{List column of Entrez IDs of all genes in each panel}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"gene_panels"

#' GENIE consortium Gene Meta Data
#'
#' Dataframe labeling all genes found in the GENIE consortium along with their corresponding
#' platform and Entrez ID.
#'
#' @format A data frame with 469 observations and 5 variables
#' \describe{
#'   \item{hugo_symbol}{Column containing all HUGO symbols genes included in IMPACT}
#'   \item{MSK_341}{Character, indicates whether each gene was included in
#'   IMPACT platform 341. Options are `included` and `not included`}
#'   \item{MSK_410}{Character, indicates whether each gene was included in
#'   IMPACT platform 410. Options are `included` and `not included`}
#'   \item{MSK_468}{Character, indicates whether each gene was included in
#'   IMPACT platform 468. Options are `included` and `not included`}
#'   \item{DFCI_1}{Character, indicates whether each gene was included in
#'   the original Dana Farber Cancer Institute sequencing panel. Options are `included` and `not included`}
#'   \item{DFCI_2}{Character, indicates whether each gene was included in
#'   the second Dana Farber Cancer Institute sequencing panel. Options are `included` and `not included`}
#'   \item{DFCI_3}{Character, indicates whether each gene was included in
#'   the third Dana Farber Cancer Institute sequencing panel. Options are `included` and `not included`}
#'   \item{DFCI_3.1}{Character, indicates whether each gene was included in
#'   the updated third Dana Farber Cancer Institute sequencing panel. Options are `included` and `not included`}
#'   \item{UHN_48}{Character, indicates whether each gene was included in
#'   the first 48 genes sequencing panel of the University Health Network. Options are `included` and `not included`}
#'   \item{UHN_50}{Character, indicates whether each gene was included in
#'   the 50 genes sequencing panel of the University Health Network. Options are `included` and `not included`}
#'   \item{VICC_1_SOLIDTUMOR}{Character, indicates whether each gene was included in
#'   the sequencing panel of the Vanderbilt-Ingram Cancer Center. Options are `included` and `not included`}
#'   \item{VICC_1_T7}{Character, indicates whether each gene was included in
#'   the sequencing panel of the Vanderbilt-Ingram Cancer Center. Options are `included` and `not included`}
#'   \item{VICC_1_T5A}{Character, indicates whether each gene was included in
#'   the sequencing panel of the Vanderbilt-Ingram Cancer Center. Options are `included` and `not included`}
#'   \item{alias}{List of EntrezIDs corresponding to the value in the hugo_symbol column and alias list.}
#'   \item{entrez_id}{List of known aliases corresponding to the value in the hugo_symbol column.}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"genie_gene_info"


#' GENIE panel names as found in synapse data
#'
#' Dataframe labeling all panels in the GENIE data and their corresponding names.
#'
#' @format A data frame with 469 observations and 5 variables
#' \describe{
#'   \item{Sequence.Assay.ID}{Column containing all panel names as found in GENIE data}
#'   \item{Panel}{Character, indicates corresponding names in gnomeR.}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"panel_names"

#' IMPACT Gene Pathways
#'
#' Dataframe labeling all genes found in IMPACT along with their corresponding
#' platform and Entrez ID.
#'
#' @format A data frame of impact genes
#' \describe{
#'   \item{pathway}{name of pathway}
#'   \item{genes}{vectors of genes in each pathways}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"pathways"



#' A segmentation file from the cbioPortal datasets
#'
#' Segmentation file provided by the processing of IMPACT sequencing using FACETS
#'
#' @format A data frame with 30240 observations with 6 variables
#' \describe{
#'  \item{ID}{Factor, IMPACT sample ID}
#'  \item{chrom}{}
#'  \item{loc.start}{}
#'  \item{loc.end}{}
#'  \item{num.mark}{}
#'  \item{seg.mean}{}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"seg"



#' Intervals sequenced in 341 panel
#'
#' List of the genes and their corresponding positions sequenced in the targeted panel of IMPACT containing 341 genes.
#' @format A data frame with 6360 lines
#' \describe{
#'     \item{V1}{Chromosome}
#'     \item{V2}{Start position}
#'     \item{V3}{End position}
#'     \item{V4}{Strand}
#'     \item{V5}{Interval name}
#'     }
#' @noRd
"ti_341"

#' Intervals sequenced in 410 panel
#'
#' List of the genes and their corresponding positions sequenced in the targeted panel of IMPACT containing 410 genes.
#'
#' @format A data frame with 7016 lines
#' \describe{
#'     \item{V1}{Chromosome}
#'     \item{V2}{Start position}
#'     \item{V3}{End position}
#'     \item{V4}{Strand}
#'     \item{V5}{Interval name}
#'     }
#' @noRd
"ti_410"


#' Intervals sequenced in 468 panel
#'
#' List of the genes and their corresponding positions sequenced in the targeted panel of IMPACT containing 468 genes.
#' @format A data frame with 7872 lines
#' \describe{
#'     \item{V1}{Chromosome}
#'     \item{V2}{Start position}
#'     \item{V3}{End position}
#'     \item{V4}{Strand}
#'     \item{V5}{Interval name}
#'     }
#' @noRd
"ti_468"

#' IMPACT Alias Tables
#'
#' Table of aliases for IMPACT panel genes
#'
#' @format A data frame with 1658 rows
#' \describe{
#'     \item{hugo_symbol}{}
#'     \item{alias}{Start }
#'     \item{entrez_id}{End }
#'     \item{alias_entrez_id}{}
#'     }
"impact_alias_table"



#' Data Frame of Column Names
#'
#' Table of accepted data names
#'
#' @format A data frame with 2333 rows
#' \describe{
#'     \item{maf_column_name}{MAD columns}
#'     \item{api_column_name}{API columns}
#'     \item{mutation_input}{Mutation columns}
#'     \item{fusion_input}{Fusion columns}
#'     \item{cna_input}{CNA columns}
#'     \item{definition}{Column Definition}
#'     \item{notes}{data notes}
#'     \item{sc_maf_column_name}{Snake case MAF name}
#'     \item{sc_api_column_name}{Snake case API column name}
#'     \item{internal_column_name}{Internally used data name}
#'     }
"names_df"



#' Consequence Map
#'
#' Dataset of recoding values for consequence mutation data
#'
#' @format A data frame with 23878 rows of CNA data from
"consequence_map"


#' An example IMPACT cBioPortal mutation data set in API format
#'
#' This set contains a random sample of 200 patients from
#' publicly available prostate cancer data from cBioPortal. The file
#' is in API format.
#'
#' @format A data frame with mutations from Abida et al. JCO Precis Oncol 2017.
#' Retrieved from cBioPortal.There are 725 observations and 29 variables.
#' \describe{
#' \item{hugoGeneSymbol}{Character w/ 324 levels,
#' Column containing all HUGO symbols genes}
#' \item{entrezGeneId}{Entrez Gene ID}
#' \item{sampleId}{MSKCC Sample ID}
#' \item{patientId}{Patient ID}
#' \item{studyId}{Indicator for Abida et al. 2017 study}
#' \item{center}{Cancer Center ID}
#' \item{mutationStatus}{Somatic or germ-line mutation status}
#' \item{variantType}{Mutation variant type}
#' \item{chr}{Chromosome mutation observed on}
#' \item{endPosition}{End Position}
#' \item{fisValue}{fisValue}
#' \item{functionalImpactScore}{functionalImpactScore}
#' \item{keyword}{keyword}
#' \item{linkMsa}{linkMsa}
#' \item{linkPdb}{linkPdb}
#' \item{linkXvar}{linkXvar}
#' \item{molecularProfileId}{molecularProfileId}
#' \item{mutationType}{mutationType}
#' \item{ncbiBuild}{ncbiBuild}
#' \item{proteinChange}{proteinChange}
#' \item{proteinPosEnd}{proteinPosEnd}
#' \item{proteinPosStart}{proteinPosStart}
#' \item{referenceAllele}{referenceAllele}
#' \item{refseqMrnaId}{refseqMrnaId}
#' \item{startPosition}{startPosition}
#' \item{uniquePatientKey}{uniquePatientKey}
#' \item{uniqueSampleKey}{uniqueSampleKey}
#' \item{validationStatus}{validationStatus}
#' \item{variantAllele}{variantAllele}}
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_mskcc_2017}
"mutations"


#' An example IMPACT cBioPortal mutation data set in API format
#'
#' This set was created from a random sample of 200 patients from
#' publicly available prostate cancer data from cBioPortal. The file
#' is in API format.
#'
#' @format A data frame with copy number alterations (CNA) from Abida et al.
#' JCO Precis Oncol 2017.Retrieved from cBioPortal.There are 475 observations
#' and 29 variables.
#' \describe{
#' \item{hugoGeneSymbol}{Character w/ 324 levels,
#' Column containing all HUGO symbols genes}
#' \item{entrezGeneId}{Entrez Gene ID}
#' \item{molecularProfileId}{Molecular Profile ID for data set}
#' \item{sampleId}{MSKCC Sample ID}
#' \item{patientId}{Patient ID}
#' \item{studyId}{Indicator for Abida et al. 2017 study}
#' \item{alteration}{Factor, Type of CNA}
#' \item{uniqueSampleKey}{character COLUMN_DESCRIPTION}
#' \item{uniquePatientKey}{character COLUMN_DESCRIPTION}}
#'
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_mskcc_2017}
"cna"


#' An example IMPACT cBioPortal mutation data set in API format
#'
#' This set was created from a random sample of 200 patients from
#' publicly available prostate cancer data from cBioPortal. The file
#' is in API format.
#'
#' @format A data frame with structural variants from Abida et al. JCO Precis Oncol 2017.
#' Retrieved from cBioPortal.There are 94 observations and 29 variables.
#' @format A data frame with 94 rows and 44 variables:
#' \describe{
#'   \item{\code{uniqueSampleKey}}{character COLUMN_DESCRIPTION}
#'   \item{\code{uniquePatientKey}}{character COLUMN_DESCRIPTION}
#'   \item{\code{molecularProfileId}}{character COLUMN_DESCRIPTION}
#'   \item{\code{sampleId}}{MSKCC Sample ID}
#'   \item{\code{patientId}}{Patient ID}
#'   \item{\code{studyId}}{Indicator for Abida et al. 2017 study}
#'   \item{\code{site1EntrezGeneId}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{site1HugoSymbol}}{Character w/ 31 levels,
#' Column containing HUGO symbols genes for first site of fusion}
#'   \item{\code{site1EnsemblTranscriptId}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site1Chromosome}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site1Position}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{site1Contig}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site1Region}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site1RegionNumber}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{site1Description}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site2EntrezGeneId}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{site2HugoSymbol}}{Character w/ 21 levels,
#' Column containing all HUGO symbols genes for second site of fusion}
#'   \item{\code{site2EnsemblTranscriptId}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site2Chromosome}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site2Position}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{site2Contig}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site2Region}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site2RegionNumber}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{site2Description}}{character COLUMN_DESCRIPTION}
#'   \item{\code{site2EffectOnFrame}}{character COLUMN_DESCRIPTION}
#'   \item{\code{ncbiBuild}}{character COLUMN_DESCRIPTION}
#'   \item{\code{dnaSupport}}{Factor, all are `yes` in this data}
#'   \item{\code{rnaSupport}}{Factor, all are `unknown` in this data}
#'   \item{\code{normalReadCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{tumorReadCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{normalVariantCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{tumorVariantCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{normalPairedEndReadCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{tumorPairedEndReadCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{normalSplitReadCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{tumorSplitReadCount}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{annotation}}{character COLUMN_DESCRIPTION}
#'   \item{\code{breakpointType}}{character COLUMN_DESCRIPTION}
#'   \item{\code{connectionType}}{character COLUMN_DESCRIPTION}
#'   \item{\code{eventInfo}}{character COLUMN_DESCRIPTION}
#'   \item{\code{variantClass}}{character COLUMN_DESCRIPTION}
#'   \item{\code{length}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{comments}}{character COLUMN_DESCRIPTION}
#'   \item{\code{svStatus}}{character COLUMN_DESCRIPTION}
#'}
#' @source \url{https://www.cbioportal.org/study/summary?id=prad_mskcc_2017}
"sv"


#' An example IMPACT cBioPortal CNA in wide format
#'
#' This set was created from a sample of 20 patients from
#' publicly available prostate cancer data from cBioPortal (`study_id = "gbc_mskcc_2022"`).
#'
#' @format A data frame with copy number alterations (CNA) retrieved from cBioPortal.
"cna_wide"

