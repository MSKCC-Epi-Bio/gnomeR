#' Public Gene Panels on cBioPortal
#'
#' Data frame of cBioPortal gene panels sourced from both public and GENIE cBioPortal instances.
#'
#' @format A nested data frame
#' \describe{
#'   \item{gene_panels}{Gene panel ID}
#'   \item{genes_in_panel}{List column of Hugo symbols of all genes in each panel}
#'   \item{entrez_ids_in_panel}{List column of Entrez IDs of all genes in each panel}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"gene_panels"


#' IMPACT Oncogenic Signaling Pathways
#'
#' Oncogenic Signaling Pathways curated from [Sanchez-Vega, F et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29625050/).
#' See [cBioPortal.org](http://www.cbioportal.org/) for more information
#'
#' @format A list of common cancer pathways and their associated alterations
#' \describe{
#'   \item{pathway}{name of pathway}
#'   \item{genes}{vector of gene alterations in each pathways}
#' }
#' @source Sanchez-Vega, F., Mina, M., Armenia, J., Chatila, W. K., Luna, A., La, K. C., Dimitriadoy, S., Liu, D. L., Kantheti, H. S., Saghafinia, S., Chakravarty, D., Daian, F., Gao, Q., Bailey, M. H., Liang, W. W., Foltz, S. M., Shmulevich, I., Ding, L., Heins, Z., Ochoa, A., … Schultz, N. (2018). Oncogenic Signaling Pathways in The Cancer Genome Atlas. Cell, 173(2), 321–337.e10. <https://doi.org/10.1016/j.cell.2018.03.035>
"pathways"


#' A segmentation file from the cbioPortal datasets
#'
#' Segmentation file provided by the processing of IMPACT sequencing using FACETS
#'
#' @format A data frame with 30240 observations with 6 variables
#' \describe{
#'  \item{ID}{Factor, IMPACT sample ID}
#'  \item{chrom}{chromosome}
#'  \item{loc.start}{start location}
#'  \item{loc.end}{end location}
#'  \item{num.mark}{number of probes or bins covered by the segment}
#'  \item{seg.mean}{segment mean value, usually in log2 scale}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"seg"


#' IMPACT Alias Tables
#'
#' Data frame of genes and their aliases for
#' IMPACT panel genes. This is used for gene name resolution functionality.
#'
#' @format A data frame with 1658 rows
#' \describe{
#'     \item{hugo_symbol}{gene Hugo Symbol}
#'     \item{alias}{Alias of Hugo Symbol in `hugo_symbol` column}
#'     \item{entrez_id}{entrez ID of gene in `hugo_symbol`}
#'     \item{alias_entrez_id}{entrez ID of `alias` gene}
#'     }
"impact_alias_table"


#' Data Frame of Column Names
#'
#' Data frame of accepted data names for standard genomic files. This serves as a
#' dictionary to help disambiguate raw column names from user entered mutation,
#' CNA or structural variant data
#'
#' @format A data frame
#' \describe{
#'     \item{maf_column_name}{data field names as they appear in common MAF file}
#'     \item{api_column_name}{data field names as they appear in common cBioPortal API retrieved files}
#'     \item{mutation_input}{does this field appear in mutation files?}
#'     \item{fusion_input}{does this field appear in mutation/sv files?}
#'     \item{cna_input}{does this field appear in CNA files?}
#'     \item{definition}{variable definition}
#'     \item{notes}{data notes}
#'     \item{sc_maf_column_name}{snake case version of `maf_column_name`}
#'     \item{sc_api_column_name}{snake case version of `api_column_name`}
#'     \item{internal_column_name}{name used for each field for all internal processing functions}
#'     }
"names_df"



#' Consequence Map
#'
#' Data frame used as a data dictionary to recode common variant classification types
#' to standardized types that can be used in oncoKB annotation.
#'
#' @format A data frame
#' \describe{
#'     \item{variant_classification}{character indicating type of mutation/variant classification as it appears in common mutation files}
#'     \item{consequence_final_coding}{final value to recode `variant_classification` column to}
#'     \item{consequence_final_coding_2}{final value to recode `variant_classification` column to}
#'     \item{consequence_final_coding_3}{final value to recode `variant_classification` column to}
#'     }
#'     @source \url{https://github.com/oncokb/oncokb-annotator/blob/a80ef0ce937c287778c36d45bf1cc8397539910c/AnnotatorCore.py#L118}
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

