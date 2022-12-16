library(cbioportalR)
library(dplyr)

#set_cbioportal_db("public")
#cna_wide <- get_cna_by_study("gbc_mskcc_2022")



cna_wide <- data.frame(
     stringsAsFactors = FALSE,
          check.names = FALSE,
          Hugo_Symbol = c("SMAD4","CCND1","MYC",
                          "FGF4","FGF3","FGF19","TAP1","ERRFI1","STK19",
                          "CRKL","SCG5","STK11","MEN1","B2M","TAP2","PMAIP1",
                          "H3-3A","H3-3B","CDC73","PIK3CA","PIK3CB",
                          "PIK3CD","PIK3CG","IGF1R","NUP93"),
  `P-0070637-T01-IM7` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,-2L),
  `P-0042589-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0026544-T01-IM6` = c(-2L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0032011-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0049337-T01-IM6` = c(0L,2L,2L,2L,2L,2L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0049685-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0022803-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0016706-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0003363-T01-IM5` = c(-2L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0068610-T01-IM7` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0068716-T01-IM7` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0022347-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0039786-T01-IM6` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L),
  `P-0065796-T01-IM7` = c(0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,0L,
                          0L,0L,0L,0L,0L)
)
usethis::use_data(cna_wide, overwrite = TRUE)
