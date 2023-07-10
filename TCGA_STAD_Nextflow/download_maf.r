library(GenomicDataCommons)
library(maftools)
library(data.table)

GenomicDataCommons::status()

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

# Get the first command-line argument
args <- commandArgs(trailingOnly = TRUE)
subdir <- args[1]

# Create a new subdirectory
if (!dir.exists(file.path(getwd(), subdir))) {
  dir.create(file.path(getwd(), subdir))
}

# Move downloaded files to the subdirectory
for (file in maffile) {
  file.copy(file, file.path(getwd(), subdir, basename(file)))
}