# context("check make oncoprint dat")
#
# test_that("make data nothing special",{
#
#   set.seed(123)
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#   bin.mut <- binmat(samples = samples,maf = mut,mut_type = "somatic_only",
#   snp_only = FALSE,include_silent = FALSE, specify_panel = FALSE)
#   gen_dat <- bin.mut[,
#   names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
#   out <- dat_oncoprint(gen_dat)
#   expect_true(is.matrix(out))
#   expect_true(all(dim(out) == c(15,200)))
#
# })
#
#
# test_that("make data while throwing warnings",{
#
#   set.seed(123)
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#   bin.mut <- binmat(samples = samples,maf = mut,mut_type = "somatic_only",
#                     snp_only = FALSE,include_silent = FALSE, specify_panel = FALSE)
#   gen_dat <- bin.mut[,
#                      names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
#   gen_dat[,1] <- sample(c("A","B"), size = nrow(gen_dat), replace = TRUE)
#   expect_warning(out <- dat_oncoprint(gen_dat))
#   expect_true(is.matrix(out))
#   expect_true(all(dim(out) == c(15,200)))
#
# })
#
#
# test_that("make data while throwing warnings",{
#
#   set.seed(123)
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#   bin.mut <- binmat(samples = samples,maf = mut,mut_type = "somatic_only",
#                     snp_only = FALSE,include_silent = FALSE, specify_panel = FALSE)
#   gen_dat <- bin.mut[,
#                      names(sort(apply(bin.mut,2, sum),decreasing = TRUE))[1:15]]
#   gen_dat[,1] <- sample(c(0:2,NA), size = nrow(gen_dat), replace = TRUE)
#   expect_warning(out <- dat_oncoprint(gen_dat))
#   expect_true(is.matrix(out))
#   expect_true(all(dim(out) == c(14,200)))
#
# })
#
#
# test_that("make data with clinical",{
#
#   set.seed(123)
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))
#   bin.mut <- binmat(samples = samples,maf = mut, fusion = fusion,
#                     cna = cna, mut_type = "somatic_only",
#                     snp_only = FALSE,include_silent = FALSE, specify_panel = FALSE)
#   gen_dat <- bin.mut[,
#                      c("TP53","CDKN2A","CDKN2A.Del","PIK3CA")]
#   gen_dat$TP53 <- rbinom(nrow(gen_dat),1,1/2)
#   gen_dat$CDKN2A.Del <- rbinom(nrow(gen_dat),1,1/2)
#   clin.patients.dat <-
#   clin.patients[match(abbreviate(rownames(gen_dat),
#   strict = TRUE, minlength = 9),
#   clin.patients$X.Patient.Identifier),] %>%
#   dplyr::rename(DMPID = X.Patient.Identifier,
#    Smoker = Smoking.History) %>%
#     dplyr::select(DMPID, Sex,Smoker) %>%
#     dplyr::filter(!is.na(DMPID)) %>%
#     dplyr::distinct(DMPID,.keep_all = TRUE) %>%
#     dplyr::mutate(Smoker = ifelse(Smoker == "Unknown",NA, as.character(Smoker)))
#   clin.patients.dat$ContFeature <- rnorm(n = nrow(clin.patients.dat))
#   gen_dat <- gen_dat[match(clin.patients.dat$DMPID,
#   abbreviate(rownames(gen_dat),strict = TRUE, minlength = 9)),]
#   clin.patients.dat <- clin.patients.dat %>%
#     tibble::rownames_to_column("to_rm") %>%
#     select(-one_of("to_rm")) %>%
#     tibble::column_to_rownames('DMPID')
#   rownames(gen_dat) <- rownames(clin.patients.dat)
#   out <- dat_oncoprint(gen_dat,clin.patients.dat)
#   expect_true(is.matrix(out))
#   expect_true(all(dim(out) == c(6,454)))
#
# })
#
#
# test_that("make data with cna non binary",{
#
#   set.seed(123)
#   samples <- as.character(unique(mut$Tumor_Sample_Barcode))[1:200]
#   bin.mut <- binmat(samples = samples,maf = mut,cna = cna, cna.binary = FALSE,
#                     mut_type = "somatic_only",
#                     snp_only = FALSE,include_silent = FALSE, specify_panel = FALSE)
#   keep <- c("TP53|PIK3CA|ALK")
#   gen_dat <- bin.mut[,grep(keep,colnames(bin.mut))]
#   expect_warning(out <- dat_oncoprint(gen_dat))
#   expect_true(is.matrix(out))
#   expect_true(all(dim(out) == c(3,200)))
#
# })
#
