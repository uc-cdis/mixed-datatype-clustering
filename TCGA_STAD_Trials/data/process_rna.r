#!/usr/bin/env Rscript
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(msigdbr))
library(BiocFileCache)
library(rtracklayer)
library(plyr)
library(edgeR)
library(org.Hs.eg.db)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
    stop("Wrong Arguments", call.=FALSE)
} else {
    raw_data_rds = args[1]
    output_path = args[2]
    cpm_threshold = args[3]
}
# remove meta feature ids
#rna_data_rds = "raw_rna_data.rds"
gencode_hyperlink = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36'
# genome name mapped to in TCGA dataset
genome_name = 'GRCh38'
# genome lengths written to RDS
gtf_rds_path = 'gene_lengths_gencode_v36.rds'
# gencode file to append to gene annotation link
gencode_file = 'gencode.v36.annotation.gtf.gz'
stad_se <- readRDS(raw_data_rds)
rownames(stad_se) <- gsub('\\.[0-9]*', '', rownames(stad_se))
stad_dge = DGEList(counts = assay(stad_se), genes = rowData(stad_se))
# filter low gene counts based on cpm threshold
prop_expressed = rowMeans(cpm(stad_dge) > 1)
#keep = prop_expressed > 0.5
keep = prop_expressed > cpm_threshold
stad_dge = stad_dge[keep, , keep.lib.sizes = FALSE]
stad_se = stad_se[keep, ]

# get gene lengths or retrieve if calculated
if (file.exists(gtf_rds_path)) {
    gene_lengths <- readRDS(gtf_rds_path)
} else {
    gencode_link = paste(gencode_hyperlink,gencode_file,sep = '/')
    bfc <- BiocFileCache()
    gencode_path <- bfcrpath(bfc, gencode_link)
    gtf = import.gff(gencode_path, format = 'gtf', genome = genome_name, feature.type = 'exon')
    grl = reduce(split(gtf, elementMetadata(gtf)$gene_id))
    gene_lengths = ldply(grl, function(x) {return(c('gene_length' = sum(width(x))))}, .id = 'ensembl_gene_id')
    genetype = unique(elementMetadata(gtf)[, c('gene_id', 'gene_type')])
    colnames(genetype)[1] = 'ensembl_gene_id'
    gene_lengths = merge(genetype, gene_lengths)
    # remove meta features similar to TCGA data
    gene_lengths$ensembl_gene_id = gsub('\\.[0-9]*', '', gene_lengths$ensembl_gene_id)
    saveRDS(gene_lengths, file = gtf_rds_path)
}
# annotate gene expression data with gene lengths
rownames(gene_lengths) = gene_lengths$ensembl_gene_id
rowData(stad_se)$gene_length = gene_lengths[rownames(stad_se), 'gene_length']
rowData(stad_se)$gene_biotype = gene_lengths[rownames(stad_se), 'gene_type']
stad_dge$genes$length = gene_lengths[rownames(stad_dge), 'gene_length']

# normalize via trimmed mean of M-values normalization method
stad_dge_tmm = calcNormFactors(stad_dge, method = 'TMM')
assay(stad_se, 'logFPKM_TMM') = rpkm(stad_dge_tmm, log = TRUE)

# add entrezid column as they are more universal
rowData(stad_se)$entrezgene = mapIds(
    org.Hs.eg.db,keys = rownames(stad_se),
    keytype = 'ENSEMBL',
    column = 'ENTREZID',
    multiVals = 'asNA'
    )
# remove duplicate entrez ids
gene_annot = rowData(stad_se)
keep = !is.na(gene_annot$entrezgene)
dup_entrez = gene_annot$entrezgene[duplicated(gene_annot$entrezgene)]
keep = keep & !gene_annot$entrezgene %in% dup_entrez
head(sort(table(gene_annot[!keep, 'gene_biotype']), decreasing = TRUE), n = 10)
stad_se = stad_se[keep,]
# save to rds file to speed things up
#saveRDS(stad_se, file = "rds_data_new_logfpkm_tmm.rds")
saveRDS(stad_se, file = output_path)
#stad_se
