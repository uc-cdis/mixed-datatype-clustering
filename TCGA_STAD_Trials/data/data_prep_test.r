#--------------------------------------------
# Kye Nichols
# Based off: Bhuva DD, Foroutan M, Xie Y, Lyu R, Cursons J, Davis MJ. Using singscore to predict mutation status in acute myeloid leukemia from transcriptomic signatures. F1000Res. 2019 Jun 3;8:776. doi: 10.12688/f1000research.19236.3. PMID: 31723419; PMCID: PMC6844140.
# original singscore publication: Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. Single sample scoring of molecular phenotypes. BMC Bioinformatics. 2018 Nov 6;19(1):404. doi: 10.1186/s12859-018-2435-4. PMID: 30400809; PMCID: PMC6219008.
# includes J's original MATH score code for the sake of simplicity
# this script calculates gene signature scores and combines with MATH and clinical data
#--------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(BiocFileCache)
library(rtracklayer)
library(plyr)
library(org.Hs.eg.db)
library(GSEABase)
library(singscore)
library(maftools)
library(data.table)
library(GenomicDataCommons)

GenomicDataCommons::status()

# get data from TCGA and process
# de_fname: name of RDS file to write expression data
process_rnaseq <- function(de_fname) {
    # gene annotation hyperlink
    gencode_hyperlink = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36'
    # genome name mapped to in TCGA dataset
    genome_name = 'GRCh38'
    # genome lengths written to RDS
    gtf_rds_path = 'gene_lengths_gencode_v36.rds'
    # gencode file to append to gene annotation link
    gencode_file = 'gencode.v36.annotation.gtf.gz'

    # check if de_fname already exists: don't overwite
    if (file.exists(de_fname)) {
        stad_se = readRDS(de_fname)
    } else {
        gdc_info = getGDCInfo()
        gdc_info
        # specify data to retrieve ftom TCGA-STAD
        # use specific sequencer
        query_rna = GDCquery(
          project = 'TCGA-STAD',
          data.category = 'Transcriptome Profiling',
          data.type = 'Gene Expression Quantification',
          workflow.type = 'STAR - Counts'
        )
        # download and prepare data to temporary directory
        rnaseq_res = getResults(query_rna)
        datapath = file.path(tempdir(), 'GDCdata')
        GDCdownload(query_rna, directory = datapath)
        stad_se = GDCprepare(query_rna, directory = datapath)
        
        # remove meta feature ids
        rownames(stad_se) <- gsub('\\.[0-9]*', '', rownames(stad_se))
        stad_dge = DGEList(counts = assay(stad_se), genes = rowData(stad_se))
        # filter low gene counts based on cpm threshold
        prop_expressed = rowMeans(cpm(stad_dge) > 1)
        keep = prop_expressed > 0.5
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
        # change to log FPKM from TMM normalized counts
        assay(stad_se, 'logFPKM_TMM') = rpkm(stad_dge_tmm, log = TRUE)
        # add entrezid column as they are more universal
        rowData(stad_se)$entrezgene = mapIds(
          org.Hs.eg.db,
          keys = rownames(stad_se),
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
        saveRDS(stad_se, file = de_fname)
    }
    stad_se
}

# calculate singsores based on stad_se
# output_fname: the output table file name
# stad_se: the dataframe containing gene expression info
# signature_paths: paths to gene signature
# signature paths can be list of one gene signature or up (index=1) and down (index=2) regulated
calc_singscore <- function(signature_name, stad_se, signature_paths) {
    # get gene signatures from list of gene signature file paths
    output_fname = paste0(signature_name, '_stad_test.csv')
    signature_paths = as.list(signature_paths)
    signature_sigs <- getBroadSets(signature_paths, membersId = 'MEMBERS_EZID')
    #signature_sigs = getBroadSets(signature_paths, membersId = 'MEMBERS_EZID')
    rownames(stad_se) = rowData(stad_se)$entrezgene
    stad_ranked = rankGenes(assay(stad_se, 'logFPKM_TMM'))
    #print(signature_sigs[[1]])
    # if gene signature list length is one calculate with no direction
    if (length(signature_sigs) == 1) {
        signature_scores = simpleScore(stad_ranked, upSet = signature_sigs[[1]], knownDirection = F)
    } else {
        signature_scores = simpleScore(stad_ranked, upSet = signature_sigs[[1]], downSet = signature_sigs[[2]])
    }
    write.csv(signature_scores, output_fname)
}

# code from J, calculating MATH score
process_mathscore <- function(math_merged_fname) {
    if (file.exists(math_merged_fname)) {
        merged_sample_summary <- readRDS(math_merged_fname)
    } else {
        cases_by_project <- cases() %>%
          facet("project.project_id") %>%
          aggregations()
        maf_files <- files() %>%
            filter(~ cases.project.project_id == "TCGA-STAD" &
                data_type == "Masked Somatic Mutation" &
                data_format == "MAF"
            ) %>%
            response_all()
        #View(maf_files)
        uid <- ids(maf_files)
        maffile <- gdcdata(uid)
        datapath = file.path(tempdir(), 'MAF_Data')

        for (file in maffile) {
            file.copy(file, file.path(datapath, basename(file)))
        }
        gz_files <- list.files(path = datapath, pattern = "*.maf.gz")
        maf_list <- list()
        # Initialize a counter for the "No non-synonymous mutations found" errors
        error_count <- 0
        # Initialize an empty data.table for storing MATH scores
        math_scores <- data.table()
        # Initialize an empty data frame for the count matrix
        count_matrix <- data.frame()
        # Process each .maf.gz file
        for (gz_file in maffile) {
          # Get the actual file path if gz_file is a symbolic link
          command <- paste("gzip -d", gz_file)
          system(command)
          # Get the .maf file name
          maf_file <- sub(".gz$", "", gz_file)
          print(maf_file)
          # Try to read the .maf file into a MAF object
          tryCatch({
            maf <- maftools::read.maf(maf = maf_file)
            # Calculate the MATH score
            math_score <- maftools::math.score(maf)
            # Merge the MATH score into the math_scores data.table
            math_scores <- rbind(math_scores, math_score)
            # Add the MAF object to the list
            maf_list[[length(maf_list) + 1]] <- maf
          },
          error = function(e) {
              # If it is, increment the error counter
              error_count <- error_count + 1
            })
          file.remove(maf_file)
        }
        # Combine the MAF objects
        merged_maf <- maftools::merge_mafs(maf_list)
        # Print the number of "No non-synonymous mutations found" errors
        print(paste("Number of 'Number of files with errors:", error_count))
        # Calculate the MATH score
        math_score_merged <- maftools::math.score(merged_maf)
        # Write the combined MAF object to a file
        datapath = file.path(tempdir(), "MAF_Merged")
        maftools::write.mafSummary(maf = merged_maf, basename = datapath)
        # Read the merged_sampleSummary.csv into a data frame
        merged_sample_summary <- read.delim(file.path(tempdir(),"MAF_Merged_sampleSummary.txt"))
        # Merge the MATH scores into the merged_sampleSummary data frame
        merged_sample_summary <- merge(merged_sample_summary, math_scores, by = "Tumor_Sample_Barcode")
        saveRDS(merged_sample_summary, file = math_merged_fname)
    }
    merged_sample_summary
}

# process clinical data and save as RDS
# pass RDS file name and write if doesn't exist
process_clinical <- function(clinical_fname) {
    # if file exists load RDS
    if (file.exists(clinical_fname)) {
        clinical <- readRDS(clinical_fname)
    } else {
        clinical <- GDCquery_clinic(project = "TCGA-STAD", type = "clinical")
        saveRDS(clinical, file = clinical_fname)
    }
    # return pulled data or saved RDS
    clinical
}

# pass as command line arguments in the future
de_fname = "rnaseq_data.rds"
math_score_rds = "MATH_SCORES.rds"
clinical_data_rds = "clinical_data.rds"

# text file containing lines with the following format:
# Gene_signature_name:gene_signature_fname_up.xml,gene_signature_fname_down.xml
# Alternatively include only one gene set regardless of up/down regulation:
# Gene_signature_name:gene_signature_fname.xml
gene_sig_config = "gene_signature_test.txt"

# process rnaseq data if doesn't exist
stad_se_in <- process_rnaseq(de_fname)

# get file names with gene signature/list and corresponding name
my_data <- readLines(gene_sig_config)
signatures <- list()
signature_dfs <- list()
for (i in my_data) {
    # gene signature name is string before colon
    tokens = strsplit(i, split = ':')
    sig_name = tokens[[1]][1]
    fname_tokens = tokens[[1]][2]
    # gene signature can be seperated by comma
    # first signature are upregulated genes in pathway
    # second signature behind comma are downregulated genes
    if (grepl(',', fname_tokens)) {
        gene_tokens <- strsplit(fname_tokens, split = ',')
        upregulated <- gene_tokens[[1]][1]
        downregulated <- gene_tokens[[1]][2]
        # gene set can be split by upregulated and downregulated genes
        gene_list <- c(upregulated, downregulated)
        
    } else {
        # no differentiation between up/down regulation
        gene_list <- c(fname_tokens)
    }
    print(gene_list)
    # call calculate singscore for every file name(s)
    calc_singscore(sig_name, stad_se_in, gene_list)
    sig_df <- read.csv(paste0(sig_name, '_stad_test.csv'), row.names=1)
    sig_df = subset(sig_df, select = c("TotalScore","TotalDispersion"))
    # include TotalDispersion in addition to total score?
    names(sig_df)[1] <- paste0(sig_name,"_Singscore")
    names(sig_df)[2] <- paste0(sig_name,"_Dispersion")
    # add dataframe to list for merging
    signature_dfs[[length(signature_dfs) + 1]] <- sig_df
    signatures[[length(signatures) + 1]] <- sig_name
}

# combine all data
singscore_merged <- data.frame(signature_dfs, stringsAsFactors=TRUE)
merged_sample_summary <- process_mathscore(math_score_rds)
clinical <- process_clinical(clinical_data_rds)

# Also J's code: 1. Define a function to extract 'bcr_patient_barcode' from 'Tumor_Sample_Barcode'
extract_bcr_patient_barcode <- function(tumor_sample_barcode) {
    # Split the string by "-"
    split_string <- strsplit(tumor_sample_barcode, "-")[[1]]
    # Concatenate the parts before the third dash
    paste(split_string[1:3], collapse = "-")
}

# Also J's code: 2. Apply the function to 'Tumor_Sample_Barcode'
merged_sample_summary$bcr_patient_barcode <- sapply(merged_sample_summary$Tumor_Sample_Barcode, extract_bcr_patient_barcode)
singscore_merged$bcr_patient_barcode <- sapply(row.names(singscore_merged), extract_bcr_patient_barcode)

# Also J's code: 3. Now you can merge 'clinical' and 'merged_sampleSummary'
# based on 'bcr_patient_barcode'
merged_df <- merge(clinical, merged_sample_summary, by = "bcr_patient_barcode")
# add singscore to merged dataframe
merged_df <- merge(merged_df, singscore_merged, by = "bcr_patient_barcode")

categorical_cols <- c("tissue_or_organ_of_origin",
                    "primary_diagnosis",
                    "prior_malignancy",
                    "ajcc_staging_system_edition",
                    "ajcc_pathologic_t",
                    "morphology",
                    "ajcc_pathologic_n",
                    "ajcc_pathologic_m",
                    "icd_10_code",
                    "site_of_resection_or_biopsy",
                    "race",
                    "ethnicity",
                    "treatments_radiation_treatment_or_therapy")

binary_cols <- c("prior_malignancy",
                "gender",
                "vital_status")

numerical_cols <- c("age_at_diagnosis",
                    "year_of_diagnosis",
                    "Frame_Shift_Del",
                    "Frame_Shift_Ins",
                    "In_Frame_Del",
                    "In_Frame_Ins",
                    "Missense_Mutation",
                    "Nonsense_Mutation",
                    "Nonstop_Mutation",
                    "Splice_Site",
                    "Translation_Start_Site",
                    "MedianAbsoluteDeviation",
                    "MATH"
                    )

# add scores from gene expression to numerical datatypes
for (sig in signatures) {
    numerical_cols[[length(numerical_cols) + 1]] <- paste0(sig, "_Singscore")
    numerical_cols[[length(numerical_cols) + 1]] <- paste0(sig, "_Dispersion")

}
all_cols <- c(categorical_cols, binary_cols, numerical_cols)
subset_df <- merged_df[, all_cols]
# for different gene signature datasets, rows are not usually removed
cleaned_df <- na.omit(subset_df)
# write to merged datafile used in run.py
write.csv(cleaned_df, "binary_numerical_merged.csv")
