#' IMPACT Gene Meta Data
#'
#' Dataframe labeling all genes found in IMPACT along with their corresponding
#' platform and Entrez ID.
#'
#' @format A data frame with 574 observations and 3 variables
#' \describe{
#'   \item{hugo_symbol}{Factor w/ 574 levels,
#'    Column containing all HUGO symbols genes included in IMPACT}
#'   \item{entrez_gene_id}{Integer, contains all Entrez IDs for genes included in IMPACT}
#'   \item{platform_added}{Character, indicates platform at which each gene was added.
#'   Options are `341`, `410`, `468`}
#'   \item{platform_341}{Character, indicates whether each gene was included in
#'   IMPACT platform 341. Options are `included` and `not included`}
#'   \item{platform_410}{Character, indicates whether each gene was included in
#'   IMPACT platform 410. Options are `included` and `not included`}
#'   \item{platform_468}{Character, indicates whether each gene was included in
#'   IMPACT platform 468. Options are `included` and `not included`}
#'   \item{pathway}{Character, indicates whether each gene belongs to one of the pathways
#'   defined in Oncogenic Signaling Pathways in The Cancer Genome Atlas`}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"impact_gene_info"


#' IMPACT Gene Meta Data
#'
#' Dataframe labeling all genes found in IMPACT along with their corresponding
#' platform and Entrez ID.
#'
#' @format A data frame of impact genes
#' \describe{
#'   \item{hugo_symbol}{Factor w/ 574 levels,
#'    Column containing all HUGO symbols genes included in IMPACT}
#'   \item{platform_341}{Character, indicates whether each gene was included in
#'   IMPACT platform 341. Options are `included` and `not included`}
#'   \item{platform_410}{Character, indicates whether each gene was included in
#'   IMPACT platform 410. Options are `included` and `not included`}
#'   \item{platform_468}{Character, indicates whether each gene was included in
#'   IMPACT platform 468. Options are `included` and `not included`}
#'   \item{pathway}{Character, indicates whether each gene belongs to one of the pathways
#'   defined in Oncogenic Signaling Pathways in The Cancer Genome Atlas`}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"impact_genes"


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


#' An example of clinical patient level information file from
#' IMPACT cbioPortal dataset.
#'
#' This set contains the publicly available data from cbioPortal patient level
#' clinical data. The file is organized with each line corresponding
#' to the information of each patient.
#'
#' @format A data frame with 497 observations on the following 6 variables
#' \describe{
#'   \item{X.Patient.Identifier}{Factor, IMPACT patient ID}
#'   \item{Sex}{Factor, Sex information: `Female`, `Male`}
#'   \item{Patient.s.Vital.Status}{Factor, indicates whether the patient is alive or deceased}
#'   \item{Smoking.History}{Factor, Indicates if person has a history of smoking}
#'   \item{Overall.Survival..Months.}{Overall survival in months}
#'   \item{Overall.Survival.Status}{Survival status}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"clin.patients"

#' An example of clinical sample level information file from IMPACT cbioPortal dataset.
#'
#' This set contains the publicly available data from cbioPortal sample level
#' clinical data. The file is organized with each line corresponding
#' to the information of each sample.
#'
#' @format A data frame with 500 observations on the following 16 variables
#' \describe{
#'   \item{X.Patient.Identifier}{Factor, IMPACT patient ID}
#'   \item{Sample.Identifier}{Factor, IMPACT sample ID}
#'   \item{Sample.Collection.Source}{Factor, was the sample collected in-house?}
#'   \item{Specimen.Preservation.Type}{Factor, indicates whether the patient is alive or deceased}
#'   \item{Specimen.Type}{Factor, Indicates if person has a history of smoking}
#'   \item{DNA.Input}{}
#'   \item{Sample.coverage}{}
#'   \item{Tumor.Purity}{}
#'   \item{Matched.Status}{Overall survival in months}
#'   \item{Sample.Type}{Survival status}
#'   \item{Primary.Tumor.Site}{Survival status}
#'   \item{Metastatic.Site}{Survival status}
#'   \item{Sample.Class}{Survival status}
#'   \item{Oncotree.Code}{Survival status}
#'   \item{Cancer.Type}{Survival status}
#'   \item{Cancer.Type.Detailed}{Survival status}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"clin.sample"


#' An example copy number alteration raw calls file from IMPACT cbioPortal dataset
#'
#' This set contains the publicly available data from cbioPortal copy number alteration
#' data. The file is organized with the first column being the gene name
#' and all other columns being samples.
#'
#' @format A data frame with 410 observations with 501 variables
#' \describe{
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"cna"

#' An example of fusions calls file from IMPACT cbioPortal dataset
#'
#' This set contains the publicly available data from cbioPortal
#' fusion data. The file is organized a classical MAF file.
#'
#' @format A data frame with 127 observations with 10 variables
#' \describe{
#'  \item{Hugo_Symbol}{Factor w/ 574 levels,
#'    Column containing all HUGO symbols genes}
#'  \item{Entrez_Gene_Id}{Entrez Gene ID (NA in this dataset)}
#'  \item{Center}{Cancer Center ID. All are `MSKCC-DMP` in this dataset}
#'  \item{Tumor_Sample_Barcode}{Factor, IMPACT sample ID}
#'  \item{Fusion}{Factor, type of fusion}
#'  \item{DNA_support}{Factor, all are `yes` in this data}
#'  \item{RNA_support}{Factor, all are `unknown` in this data}
#'  \item{Method}{All are NA in this data}
#'  \item{Frame}{Factor, `in frame`, `out of frame` or `unknown`}
#'  \item{Comments}{}
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"fusion"

#' An example maf file from IMPACT cbioPortal dataset
#'
#' This set contains the publicly available data from cbioPortal mutation data.
#' The file is organized a classical MAF file.
#'
#' @format A data frame with 457 observations with 45 variables
#' \describe{
#' }
#' @source \url{https://cbioportal.mskcc.org/}
"mut"


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

#' List of curated genes for IMPACT oncoKB annotation.
#'
#' We include the curated genes for the IMPACT sequencing for the oncoKB annotation pipeline.
#'
#' @format A data frame with 417 observations on the following 2 variables.
#' \describe{
#'  \item{id}{}
#'  \item{gene}{}
#' }
"curated_genes"


#' List of curated genes for IMPACT
#'
#' We include the curated genes for the IMPACT sequencing from 3 different
#' platforms. This set of organized as a list of three different character
#' vectors corresponding to each of the 341, 410 and 468 gene platforms.
#'
#' @format A list of 3 items.
#' \describe{
#'  \item{g341}{}
#'  \item{g410}{}
#'  \item{g468}{}
#' }
"g.impact"



#' A vector of 19441 hugo symbols in TCGA
#'
#'A vector of 19441 hugo symbols in TCGA
#'
#' @format A data frame with 19441 genes
#' \describe{
#'     }
"tcga_genes"


#' Data frame of all TCGA sample ids and their corresponding cancer type
#'
#' Data frame of all TCGA sample ids and their corresponding cancer type
#'
#' @format A data frame with 10275 samples
#' \describe{
#'     \item{patient_id}{Sample ID}
#'     \item{Cancer_Code}{Cancer site following TCGA abbreviations}
#'     }
"tcga_samples"


#' Intervals sequenced in 341 panel
#'
#' List of the genes and their corresponding positions sequenced in the targeted panel of IMPACT containing 341 genes.
#'
#' @format A data frame with 6360 lines
#' \describe{
#'     \item{V1}{Chromosome}
#'     \item{V2}{Start position}
#'     \item{V3}{End position}
#'     \item{V4}{Strand}
#'     \item{V5}{Interval name}
#'     }
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
"ti_410"


#' Intervals sequenced in 468 panel
#'
#' List of the genes and their corresponding positions sequenced in the targeted panel of IMPACT containing 468 genes.
#'
#' @format A data frame with 7872 lines
#' \describe{
#'     \item{V1}{Chromosome}
#'     \item{V2}{Start position}
#'     \item{V3}{End position}
#'     \item{V4}{Strand}
#'     \item{V5}{Interval name}
#'     }
"ti_468"
