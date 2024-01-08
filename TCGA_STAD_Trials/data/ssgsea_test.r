library("pathwayPCA")
library(SummarizedExperiment)
library(GSEABase)
library(corto)


ssgsea_calc <- function(de_path, gmt_path) {
    # read stad_se from singscore script
    stad_se <- readRDS(de_path)
    # get TPM values from RSEM output
    ssgsea_tpm <- assay(stad_se, "tpm_unstrand")
    ssgsea_input <- data.matrix(ssgsea_tpm)
    gmt_list <- read_gmt(gmt_path, description = TRUE)
    ssgsea_mat <-ssgsea(ssgsea_input,gmt_list$pathways)
    ssgsea_df <- as.data.frame(ssgsea_mat)
    # reapply labels to output matrix
    rownames(ssgsea_df) = gmt_list$TERMS
    colnames(ssgsea_df) = colnames(ssgsea_tpm)
    ssgsea_df
}


ssgsea_output <- ssgsea_calc("stad_de.rds", "kegg_hsa.gmt")
write.csv(ssgsea_output, "ssgsea_test_output.csv")
