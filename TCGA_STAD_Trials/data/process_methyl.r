# packages for gene set retr
library(plyr)
#library(msigdbr)

# packages for TCGA retr
library(TCGAbiolinks)
library(SummarizedExperiment)
library(BiocFileCache)
#library(maftools)
library(data.table)
#library(rtracklayer)
library(GenomicDataCommons)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
    stop("Wrong Arguments", call.=FALSE)
} else {
    methyl_rds = args[1]
    output_csv = args[2]
}

# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html
# adapted from NBIS workshop documentation
process_methyl <- function(methyl_rds_path) {
        #url =
"https://github.com/hamidghaedi/Methylation_Analysis/blob/master/cross_reactive_probe.chen2013.csv"
        #download.file(url, "cross_reactive_probe.chen2013.csv", mode = "wb")
        met <- readRDS(methyl_rds_path)
        # get probe annotation
        ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        ## remove NA probes
        probe.na <- rowSums(is.na(met))
        table(probe.na == 0)
        probe <- probe.na[probe.na == 0]
        met <- met[row.names(met) %in% names(probe), ]
        keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
        table(keep)
        met <- met[keep, ]
        rm(keep)
        ## remove SNPs overlapped probe
        table (is.na(ann450k$Probe_rs))
        # probes without snp
        no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]
        snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
        #snps with maf <= 0.05
        snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]
        # filter met
        met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]
        #remove no-further needed dataset
        rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)
        #crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
        
        #probe_file_path <- "my_methyl_probes.txt"
        #crs.reac <- read.csv(probe_file_path)
        #crs.reac <- crs.reac$TargetID[-1]
        #print(crs.reac)
        # filter met
        #met <- met[ -which(row.names(met) %in% crs.reac), ]
        #bval <- met
        ## converting beta values to m_values
        ## m = log2(beta/1-beta)
        mval <- t(apply(met, 1, function(x) log2(x/(1-x))))
        #print(mval)
        #print(met)
        #saveRDS(mval, file = methyl_rds_path, compress = FALSE)
        #methyl_rds_file <- mval
        write.csv(mval, output_csv)
}

process_methyl(methyl_rds)
