library(GenomicDataCommons)
library(maftools)
library(data.table)

GenomicDataCommons::status()

# Read in the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]
subdir <- args[2]

cases_by_project <- cases() %>%
  facet("project.project_id") %>%
  aggregations()

maf_files <- files() %>%
    filter(~ cases.project.project_id == project_name &
        data_type == "Masked Somatic Mutation" &
        data_format == "MAF"
    ) %>%
    response_all()

uid <- ids(maf_files)

maffile <- gdcdata(uid)

# Create a new subdirectory
if (!dir.exists(file.path(getwd(), subdir))) {
  dir.create(file.path(getwd(), subdir))
}

# Move downloaded files to the subdirectory
for (file in maffile) {
  file.copy(file, file.path(getwd(), subdir, basename(file)))
}