library(maftools)
library(data.table)
library(TCGAbiolinks)

# Get the first command-line argument
args <- commandArgs(trailingOnly = TRUE)
subdir <- args[1]

# Get list of .maf.gz files in the subdirectory
gz_files <- list.files(path = subdir, pattern = "*.maf.gz")

# Initialize an empty list to store the MAF objects
maf_list <- list()

# Initialize a counter for the "No non-synonymous mutations found" errors
error_count <- 0

# Initialize an empty data.table for storing MATH scores
math_scores <- data.table()

# Initialize an empty data frame for the count matrix
count_matrix <- data.frame()

# Process each .maf.gz file
for (gz_file in gz_files) {
  # Get the actual file path if gz_file is a symbolic link
  if (file.info(file.path(subdir, gz_file))$isdir == FALSE) {
    command <- paste("gzip -d $(readlink -f", file.path(subdir, gz_file), ")")
    print(command)
  } else {
    # Decompress the .maf.gz file
    command <- paste("gzip -d", file.path(subdir, gz_file))
  }
  system(command)

  # Get the .maf file name
  maf_file <- sub(".gz$", "", gz_file)

  print(maf_file)

  # Try to read the .maf file into a MAF object
  tryCatch({
    maf <- maftools::read.maf(maf = file.path(subdir, maf_file))

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

#   # If the file is a symbolic link, resolve to actual file path
#   if (file.info(file.path(subdir, maf_file))$isdir == FALSE) {
#     actual_file_path <- system(paste("readlink -f", file.path(subdir, maf_file)), intern = TRUE)
#     file.remove(actual_file_path)
#   } else {
#     file.remove(file.path(subdir, maf_file))
#   }
  file.remove(file.path(subdir, maf_file))

}

# Combine the MAF objects
merged_maf <- maftools::merge_mafs(maf_list)

# Print the number of "No non-synonymous mutations found" errors
print(paste("Number of 'Number of files with errors:", error_count))

# Calculate the MATH score
math_score_merged <- maftools::math.score(merged_maf)

# Write the combined MAF object to a file
maftools::write.mafSummary(maf = merged_maf, basename = "maf_files/merged")

# Get list of .txt files in the subdirectory
txt_files <- list.files(path = subdir, pattern = "\\.txt$")

# Process each .txt file
for (txt_file in txt_files) {
  # Read the .txt file into a data frame
  data <- read.delim(file.path(subdir, txt_file))

  # Write the data frame to a .csv file
  write.csv(data, file = file.path(subdir,
                sub("\\.txt$", ".csv", txt_file)), row.names = FALSE)
}

# Read the merged_sampleSummary.csv into a data frame
merged_sample_summary <- read.csv("maf_files/merged_sampleSummary.csv")

# Merge the MATH scores into the merged_sampleSummary data frame
merged_sample_summary <- merge(merged_sample_summary, math_scores,
                                            by = "Tumor_Sample_Barcode")

clinical <- GDCquery_clinic(project = "TCGA-STAD", type = "clinical")

# 1. Define a function to extract 'bcr_patient_barcode' from 'Tumor_Sample_Barcode'
extract_bcr_patient_barcode <- function(tumor_sample_barcode) {
    # Split the string by "-"
    split_string <- strsplit(tumor_sample_barcode, "-")[[1]]
    # Concatenate the parts before the third dash
    paste(split_string[1:3], collapse = "-")
}

# 2. Apply the function to 'Tumor_Sample_Barcode' 
# column in 'merged_sampleSummary'
merged_sample_summary$bcr_patient_barcode <- sapply(merged_sample_summary$Tumor_Sample_Barcode, extract_bcr_patient_barcode)

# 3. Now you can merge 'clinical' and 'merged_sampleSummary'
# based on 'bcr_patient_barcode'
merged_df <- merge(clinical, merged_sample_summary, by = "bcr_patient_barcode")

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
                    "MATH")

all_cols <- c(categorical_cols, binary_cols, numerical_cols)
subset_df <- merged_df[, all_cols]

# 12 rows are removed here mostly because they don't have an age at diagnosis 
# we may consider alternative ways to get around this missing data problem
cleaned_df <- na.omit(subset_df)

# Write the merged_sampleSummary data frame with the MATH scores
# back to a .csv file
write.csv(cleaned_df, "maf_files/binary_numerical_merged.csv")
# maf_files/merged_sampleSummary_with_MATH.csv

# # Generate the count matrix of mutations
# mut_count_matrix <- maftools::mutCountMatrix(merged_maf)

# # Reformat the count matrix to have samples as rows,
# # genes as columns, and binary values
# # reformatted_count_matrix = copy(as.data.frame(t(count_matrix)))

# # print(reformatted_count_matrix)
# # print(class(reformatted_count_matrix))

# write.csv(as.data.frame(t(mut_count_matrix)),
#                 "maf_files/reformatted_count_matrix.csv",  row.names = TRUE)

# reformatted_count_matrix <- read.csv("maf_files/reformatted_count_matrix.csv")
# tumor_barcodes <- reformatted_count_matrix$X

# reformatted_count_matrix <- as.data.frame(
#             lapply(reformatted_count_matrix, function(x) ifelse(x > 0, 1, 0)))

# reformatted_count_matrix$Tumor_Sample_Barcode <- tumor_barcodes

# reformatted_count_matrix <- subset(reformatted_count_matrix, select = -X)

# merged_df <- merge(reformatted_count_matrix,
#                     merged_sample_summary,
#                     by = "Tumor_Sample_Barcode", all = TRUE)

# # Drop rows with NaN values
# merged_df <- merged_df[complete.cases(merged_df), ]

# # Write the merged_sampleSummary data frame with the MATH scores
# # back to a .csv file
# write.csv(merged_df, "maf_files/binary_numerical_merged.csv")