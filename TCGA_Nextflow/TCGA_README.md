# TCGA Nextflow Workflow

This code repository provides a comprehensive workflow for analyzing clinical data and mutation annotation files (MAF files) from TCGA projects. It has been specifically tested on TCGA-PRAD and TCGA-BRCA.

## Table of Contents
- [Requirements](#requirements)
- [Usage](#usage)
- [Configuration Options](#configuration-options)
- [Scripts and Files Description](#scripts-and-files-description)

## Requirements

Have [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) in your path and have [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.

## Usage

### Running the Nextflow Workflow
Run the main workflow with the following command:

```bash
nextflow run main.nf
```

You can specify various parameters by adding them to the command line. For example:

```bash
nextflow run main.nf --project_name TCGA-STAD --categorical_cols "column1,column2" --binary_cols "column3" --numerical_cols "column4"
```

## Configuration Options

The `nextflow.config` file contains options that can be configured for the workflow. Below are the available parameters:

- `download_maf_script`: Name of the script for downloading MAF files.
- `maftools_script`: Name of the script for collecting clinical data, processing maf files, etc.
- `py_script`: Name of the script to perform clustering.
- `project_name`: Name of the TCGA project to scrape data from (e.g., "TCGA-STAD").
- `categorical_cols`: Comma-separated string listing categorical columns for the dataframe.
- `binary_cols`: Comma-separated string listing binary columns for the dataframe.
- `numerical_cols`: Comma-separated string listing numerical columns for the dataframe.
- `include_gene_level`: Indicates whether to include gene-level mutation data ("Yes" or "No").
- `dataframe_size`: Maximum number of rows * columns allowed in the dataframe for clustering.
- `subdirectory_name`: Subdirectory to put downloaded files in.
- `hyperparams`: List of hyperparameters for clustering.
- `distance_types`: List of distance types for clustering.
- `clusters`: List of numbers specifying cluster sizes.

## Scripts and Files Description

- `download_maf.r`: Downloads MAF files from specific projects.
- `maftools_processing.r`: Processes MAF files using maftools.
- `main.nf`: Main Nextflow workflow for data processing.
- `nextflow.config`: Configuration for the Nextflow workflow.
- `mixed_datatype_clustering.py`: Clustering functions for mixed data types.
- `process_data.py`: Processes expression data and visualizations.
- `python_environment.yml`: Python environment configuration.
- `r_environment.yml`: R environment configuration.