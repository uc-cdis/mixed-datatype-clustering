library(plyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(BiocFileCache)
library(data.table)
library(GenomicDataCommons)
library(SummarizedExperiment)



# get data from TCGA and process
# de_fname: name of RDS file to write expression data
download_rnaseq <- function(output_fname) {
    # if file exists load RDS
    if (file.exists(output_fname)) {
        stad_se <- readRDS(output_fname)
    } else {
        # gene annotation hyperlink
        gencode_hyperlink =
    'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36'
        # genome name mapped to in TCGA dataset
        genome_name = 'GRCh38'
        # genome lengths written to RDS
        gtf_rds_path = 'gene_lengths_gencode_v36.rds'
        # gencode file to append to gene annotation link
        gencode_file = 'gencode.v36.annotation.gtf.gz'
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
        saveRDS(stad_se, rna_out_path)
    }
    stad_se
}
    

# process clinical data and save as RDS/CSV
# pass RDS file path and write if doesn't exist
download_clinical <- function(clinical_file_path) {
    # if file exists load RDS
    if (file.exists(clinical_file_path)) {
        clinical <- readRDS(clinical_file_path)
    } else {
        clinical <- GDCquery_clinic(project = "TCGA-STAD", type = "clinical")
        #saveRDS(clinical, file = clinical_file_path)
        write.csv(clinical, clinical_file_path)
    }
    # return pulled data or saved RDS
    clinical
}

download_methyl <- function(methyl_rds_path) {
    if (file.exists(methyl_rds_path)) {
        met <- readRDS(methyl_rds_path)
    } else {
        query_met <- GDCquery(project= "TCGA-STAD",
                                   data.category = "DNA Methylation",
                                   platform = "Illumina Human Methylation 450",
                                   data.type = "Methylation Beta Value")
#                                   legacy = TRUE)
        GDCdownload(query_met)
        data.met <- GDCprepare(query_met)
        met <- as.data.frame(SummarizedExperiment::assay(data.met))
        saveRDS(met, file = methyl_rds_path, compress = FALSE)
    }
    met
}

# code from J, calculating MATH score
process_mathscore <- function(math_merged_path) {
    if (file.exists(math_merged_path)) {
        merged_sample_summary <- readRDS(math_merged_path)
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
        #saveRDS(merged_sample_summary, file = math_merged_path)
        write.csv(merged_sample_summary, ath_merged_path)
    }
    merged_sample_summary
}


# https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/subtypes.html
# save labeling subtypes
get_stad_subtypes <- function(subtype_path) {
    if (file.exists(subtype_path)) {
        subtype_data <- readRDS(subtype_path)
    } else {
        subtypes <- PanCancerAtlas_subtypes()
        subtype_data <- DT::datatable(
            data = subtypes,
            filter = 'top',
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
            rownames = FALSE
            )
        #saveRDS(subtype_data, file = subtype_path)
        subtypes_data = subtypes_data[[1]]$data
        write.csv(subtypes_data, subtype_path)
    }
    subtype_data
}


args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
    stop("Wrong Arguments", call.=FALSE)
} else {
    rna_out_path = args[1]
    download_rnaseq(rna_out_path)

    clinical_out_path = args[2]
    download_clinical(clinical_out_path)

    methyl_out_path = args[3]
    download_methyl(methyl_out_path)
    
    math_out_path = args[4]
    process_mathscore(math_out_path)
    
    subtype_out_path = args[5]
    get_stad_subtypes(subtype_out_path)

}

