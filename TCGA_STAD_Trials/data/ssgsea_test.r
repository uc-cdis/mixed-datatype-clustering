#!/usr/bin/env Rscript
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
#suppressMessages(library(circlize))
#suppressMessages(library(corto))
suppressMessages(library(msigdbr))
#source("ssgsea_rpubs.r")

# ssGsea Rpubs implimentation: https://rpubs.com/pranali018/SSGSEA
# exact function included here as ssgsea_rpub
ssgsea_rpub = function(X, gene_sets, alpha = 0.25, scale = T, norm = T, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)
        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos
            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
            step_cdf_diff = step_cdf_pos - step_cdf_neg
            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes
            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })
    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))
    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

args = commandArgs(trailingOnly=TRUE)
print(args)
print(length(args))
if (length(args)!=7) {
    stop("Wrong Arguments", call.=FALSE)
} else {
    gs_cat_fname = args[1]
    msigdb_cats = strsplit(args[2], ",")
    scaled_tf = args[3]
    norm_tf = args[4]
    rna_seq_rds = args[5]
    rna_metric = args[6]
    output_csv = args[7]
}

my_gene_sets  <- readLines(gs_cat_fname)
geneset_list <- c()
for (msigdb_cat in msigdb_cats) {
    gs = msigdbr(species = "human", category = msigdb_cats)
    geneset_list <- c(geneset_list, gs)
}
genesets = dplyr::bind_rows(geneset_list)
genesets = genesets %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
gs_data <- genesets[genesets$gs_name %in% my_gene_sets,]
gs_data <-split(x = gs_data$ensembl_gene, f = gs_data$gs_name)
stad_se = readRDS(rna_seq_rds)

rna_data <- assay(stad_se, rna_metric)
log2_mat <- data.matrix(log2(rna_data)+1)
input_df <- as.data.frame(log2_mat)

#ssgsea_mat <-ssgsea(input_df, gs_data, scale = scaling, minsize = minsize)
if (scaled_tf == "T") {
    scaled_bool = T
} else {
    scaled_bool = F
}
if (norm_tf == "T") {
    norm_bool = T
} else {
    norm_bool = F
}
ssgsea_mat <-ssgsea_rpub(data.matrix(input_df), gs_data, alpha = 0.25, scale = scaled_bool, norm = norm_bool, single = T)
ssgsea_df <- as.data.frame(ssgsea_mat)
colnames(ssgsea_df) = colnames(rna_data)
write.csv(ssgsea_df, output_csv)
