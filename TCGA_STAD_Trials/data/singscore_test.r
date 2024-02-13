#--------------------------------------------
# Kye Nichols
# Based off: Bhuva DD, Foroutan M, Xie Y, Lyu R, Cursons J, Davis MJ. Using singscore to predict mutation status in acute myeloid leukemia from transcriptomic signatures. F1000Res. 2019 Jun 3;8:776. doi: 10.12688/f1000research.19236.3. PMID: 31723419; PMCID: PMC6844140.
# original singscore publication: Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. Single sample scoring of molecular phenotypes. BMC Bioinformatics. 2018 Nov 6;19(1):404. doi: 10.1186/s12859-018-2435-4. PMID: 30400809; PMCID: PMC6219008.
# includes J's original MATH score code for the sake of simplicity
# this script calculates gene signature scores and combines with MATH and clinical data

#!/usr/bin/env Rscript
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(msigdbr))
suppressMessages(library(singscore))
suppressMessages(library(msigdbr))
suppressMessages(library(BiocFileCache))
suppressMessages(library(rtracklayer))
suppressMessages(library(plyr))
suppressMessages(library("stringr"))
library(edgeR)
library(org.Hs.eg.db)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
    stop("Wrong Arguments", call.=FALSE)
} else {
    rna_data_rds = args[1]
    gs_cat_fname = args[2]
    msigdb_cats = strsplit(args[3], ",")
    inc_direction = args[4]
    rna_metric = args[5]
    output_fname = args[6]
}

up_str = "_UP"
down_str = "_DN"
my_gene_sets  <- readLines(gs_cat_fname)
geneset_list <- c()
for (msigdb_cat in msigdb_cats) {
    gs = msigdbr(species = "human", category = msigdb_cats)
    geneset_list <- c(geneset_list, gs)
}
genesets = dplyr::bind_rows(geneset_list)
genesets = genesets %>% dplyr::distinct(gs_name, entrez_gene) %>% as.data.frame()
gs_data <- genesets[genesets$gs_name %in% my_gene_sets,]
gs_data <-split(x = gs_data$entrez_gene, f = gs_data$gs_name)

calc_singscore <- function(signature, stad_se, sig_up, sig_down=NULL) {
    rownames(stad_se) = rowData(stad_se)$entrezgene
    # stad_ranked = rankGenes(assay(stad_se, 'logFPKM_TMM'))
    stad_ranked = rankGenes(assay(stad_se, rna_metric))
    if (is.null(sig_down)) {
        #print(signature)
        signature_scores = simpleScore(stad_ranked, upSet = sig_up, knownDirection = FALSE)
    } else {
        if (inc_direction=="T") {
            signature_scores = simpleScore(stad_ranked, upSet = sig_up, downSet = sig_down)
        } else {
            signature_scores = simpleScore(stad_ranked, upSet = sig_up, knownDirection = FALSE)
        }
    }
    colnames(signature_scores) <- paste(colnames(signature_scores),signature,sep="_")
    signature_scores
}

singscore_dataframes = NULL
counter = 0 #rna_data_rds = "rds_data_new_logfpkm_tmm.rds"
stad_se <- readRDS(rna_data_rds)

for (gs in as.list(names(gs_data))) {
    if (endsWith(gs, up_str) | endsWith(gs, down_str)) {
        if (endsWith(gs, up_str)) {
            down_sig = str_replace_all(gs, up_str, down_str)
            ss <-calc_singscore(gs, stad_se, gs_data[[gs]], sig_down=gs_data[[down_sig]])
            if (counter == 0) {
                singscore_dataframes = ss
            } else {
                singscore_dataframes <- dplyr::bind_cols(singscore_dataframes, ss)
            }
            counter = counter+ 1
        }
    } else {
        ss <-calc_singscore(gs, stad_se, gs_data[[gs]])
        if (counter == 0) {
            singscore_dataframes = ss
        } else {
            singscore_dataframes <- dplyr::bind_cols(singscore_dataframes, ss)
        }
        counter = counter+ 1
    }
}
write.csv(singscore_dataframes, output_fname)
