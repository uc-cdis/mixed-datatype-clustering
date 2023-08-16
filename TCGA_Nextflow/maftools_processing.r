library(maftools)
library(data.table)
library(TCGAbiolinks)

# Get the first command-line argument
args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]
subdir <- args[2]
dataframe_size <- as.numeric(args[3])
categorical_cols <- unlist(strsplit(args[4], split = ",")[[1]])
binary_cols <- unlist(strsplit(args[5], split = ",")[[1]])
numerical_cols <- unlist(strsplit(args[6], split = ",")[[1]])
include_gene_level <- args[7]

all_cols <- c(categorical_cols,
            binary_cols,
            numerical_cols,
            c("Tumor_Sample_Barcode"))

print(all_cols)

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

  file.remove(file.path(subdir, maf_file))

}

# Combine the MAF objects
merged_maf <- maftools::merge_mafs(maf_list)

# Print the number of "No non-synonymous mutations found" errors
print(paste("Number of 'Number of files with errors:", error_count))

# Calculate the MATH score
math_score_merged <- maftools::math.score(merged_maf)

# Write the combined MAF object to a file
maftools::write.mafSummary(
    maf = merged_maf, basename = paste(subdir, "/merged", sep=""))

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
merged_sample_summary <- read.csv(
    paste(subdir, "/merged_sampleSummary.csv", sep = ""))

# Merge the MATH scores into the merged_sample_summary data frame
merged_sample_summary <- merge(merged_sample_summary, math_scores,
                                            by = "Tumor_Sample_Barcode")

clinical <- GDCquery_clinic(project = project_name, type = "clinical")

# 1. Define a function to extract 'bcr_patient_barcode' from 'Tumor_Sample_Barcode'
extract_bcr_patient_barcode <- function(tumor_sample_barcode) {
    # Split the string by "-"
    split_string <- strsplit(tumor_sample_barcode, "-")[[1]]
    # Concatenate the parts before the third dash
    paste(split_string[1:3], collapse = "-")
}

# 2. Apply the function to 'Tumor_Sample_Barcode' column in 'merged_sample_summary'
merged_sample_summary$bcr_patient_barcode <- sapply(merged_sample_summary$Tumor_Sample_Barcode, extract_bcr_patient_barcode)

# 3. Now you can merge 'clinical' and 'merged_sample_summary' based on 'bcr_patient_barcode'
non_bin_merged_df <- merge(clinical, merged_sample_summary, by = "bcr_patient_barcode")

missing_cols <- all_cols[!(all_cols %in% names(non_bin_merged_df))]

if (length(missing_cols) > 0) {
  stop("These columns are not in the dataframe:", paste(missing_cols, collapse = ", "))
}

subset_df <- non_bin_merged_df[, all_cols]

# Generate the count matrix of mutations
mut_count_matrix <- maftools::mutCountMatrix(merged_maf)

write.csv(as.data.frame(t(mut_count_matrix)),
        paste(subdir, "/reformatted_count_matrix.csv", sep = ""),
        row.names = TRUE)

reformatted_count_matrix <- read.csv(
    paste(subdir, "/reformatted_count_matrix.csv", sep = ""))

tumor_barcodes <- reformatted_count_matrix$X

reformatted_count_matrix <- as.data.frame(
            lapply(reformatted_count_matrix, function(x) ifelse(x > 0, 1, 0)))

reformatted_count_matrix$Tumor_Sample_Barcode <- tumor_barcodes

reformatted_count_matrix <- subset(reformatted_count_matrix, select = -X)

if (include_gene_level == "Yes") {
    merged_df <- merge(reformatted_count_matrix,
                        subset_df,
                        by = "Tumor_Sample_Barcode", all = TRUE)
} else {
    merged_df <- subset_df
}

# Find the number of NaN's in each column
na_count <- na_count <- colSums(is.na(merged_df))

# Drop rows with NaN values
merged_df <- merged_df[complete.cases(merged_df), ]

# Check if there are no rows left
if (nrow(merged_df) == 0) {
  # List the columns with the most NaN's along with how many NaN's they have
  sorted_na_count <- sort(na_count, decreasing = TRUE)
  cat("Error: No rows left after removing NaN values. Columns with the most NaN's:\n")
  print(sorted_na_count)
  stop()
}

# Repeat until the number of rows * columns is less than 1.5 million
# This makes sure that the clustering algorithm will proceed sufficiently fast
while (nrow(merged_df) * ncol(merged_df) > dataframe_size) {
  # Take 90% of the original rows
  merged_df <- merged_df[sample(nrow(merged_df), nrow(merged_df) * 0.9),]

  # Remove all columns which contain duplicate values or where all values are the same
  merged_df <- merged_df[, sapply(merged_df, function(col) length(unique(col)) > 1)]

  # Print the current number of rows * columns
  cat("Number of rows * columns:", nrow(merged_df) * ncol(merged_df), "\n")
}

# Write the merged_sample_summary data frame with the MATH scores
# back to a .csv file
write.csv(merged_df, paste(subdir, "/binary_numerical_merged.csv", sep = ""))