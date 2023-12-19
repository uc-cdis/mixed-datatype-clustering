import gseapy
import os

# import sys
import numpy as np
import pandas as pd
import csv
import shutil

# run GSEA. ex: gseapy.gsea(data='expression.txt', gene_sets='gene_sets.gmt', cls='test.cls', outdir='test')
#      run prerank. ex: gseapy.prerank(rnk='gsea_data.rnk', gene_sets='gene_sets.gmt', outdir='test')
# run GSVA: ex: gseapy.gsva(data="expression.txt", gene_sets= "gene_sets.gmt", outdir='test')

# Resource: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA
# ssGEA operates independent of phenotypic labeling which is suitable for our analysis
# run ssGSEA. ex: gseapy.ssgsea(data="expression.txt", gene_sets= "gene_sets.gmt", outdir='test')

# gene_counts_fname is file name with gene counts (TCGA open files)
# gene_set_fname is a Gene Matrix Transposed file (*.gmt)
# myoutput_dir is the name of the output directory

# Single-sample GSEA
def ssgsea_analysis(gene_set_fname, indir, outdir):
    # file extension of input
    ext = "tsv"
    qc_dict = {}
    # summary column: if unmapped or multimapped features are not optimal we need to get access to
    # fastq/bam files and map in a way that's more optimal for our analysis
    row_omit = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]
    cols_keep = ["gene_name", "tpm_unstranded"]
    # cols_keep = ["tpm_unstranded"]

    # set paths of source and output data and gene set(s)
    wk_dir = os.getcwd()
    data_path = os.path.join(wk_dir, indir)
    out_path = os.path.join(wk_dir, outdir)
    gmt_path = os.path.join(wk_dir, gene_set_fname)

    # create output folder, replace if exists
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    else:
        shutil.rmtree(out_path)
        os.mkdir(out_path)

    temp_folder = "temp_ssgsea"
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)
    else:
        shutil.rmtree(temp_folder)
        os.mkdir(temp_folder)

    # gene count paths as a list
    gcount_paths = [
        os.path.join(data_path, i)
        for i in os.listdir(data_path)
        if i.endswith(".%s" % ext)
    ]
    # traverse gene counts
    for gcount_path in gcount_paths:
        # replace symbolic link with absolute path
        if os.path.realpath(gcount_path):
            gcount_path = os.path.realpath(gcount_path)
        if os.path.realpath(gmt_path):
            gmt_path = os.path.realpath(gmt_path)

        # rows 2,3,4,5 are qc information: omit by row_omit cols
        df = pd.read_csv(
            gcount_path, delimiter="\t", skiprows=[0, 2, 3, 4, 5], index_col=0
        )
        # alternatively use df.drop(row_omit)
        # gene type can be transcribed_unitary_pseudogene, processed pseudo gene, etc.
        # for the sake of simplicity include only protein coding
        protein_coding = df[df["gene_type"] == "protein_coding"]

        # protein_coding = protein_coding.reset_index().set_index("gene_name")
        protein_coding = protein_coding[cols_keep]
        # protein_coding.index = protein_coding.index.str.replace(r'[^.]+$', '', regex=True)
        # protein_coding["gene_name"] = protein_coding["gene_name"].str.extract(r'(\w+)', expand=True)
        # protein_coding.index = protein_coding.index.str.extract(r'(\w+)', expand=True)
        protein_coding.index = protein_coding.index.str.split(".").str[0]
        # if low number of duplicates, ssgsea will average them together

        # print(protein_coding[protein_coding["gene_name"].duplicated() == True])
        # 86/20,000 genes might be negligable
        # print(len(protein_coding[protein_coding.index.duplicated()]))

        # create seperate folder for ssgsea output with patient barcode
        sample_name = os.path.basename(gcount_path).replace(".%s" % ext, "")
        input_path = os.path.join(temp_folder, sample_name + (".txt"))

        protein_coding.to_csv(input_path, sep="\t")
        sample_folder = os.path.join(out_path, sample_name)
        if not os.path.isdir(sample_folder):
            os.mkdir(sample_folder)
        else:
            shutil.rmtree(sample_folder)
            os.mkdir(sample_folder)

        # print("Processing %s" % sample_name)
        # min_s = 1 # default is 15
        # max_s = 5000 # max is 500
        # ssgea: single sample gea -> normalized
        gseapy.ssgsea(
            data=input_path,
            sample_norm_method="rank",
            gene_sets=gmt_path,
            #            min_size = min_s,
            #            max_size = max_s,
            outdir=sample_folder,
        )
    shutil.rmtree(temp_folder)


# merge single sample analyses
def ssgsea_merge(ssgea_merged_fname, output_dir):
    wk_dir = os.getcwd()
    ssgea_merged_path = os.path.join(wk_dir, ssgea_merged_fname)
    gsea_report_fname = "gseapy.gene_set.ssgsea.report.csv"
    subdirectories = [dir for dir in os.listdir(output_dir) if "" in dir]
    sample_dict = {}
    for ssample in subdirectories:
        ssample_path = os.path.join(output_dir, ssample)
        report_path = os.path.join(ssample_path, gsea_report_fname)
        sample_dict[
            ssample.replace(".rna_seq.augmented_star_gene_counts", "")
        ] = report_path
    # create dictioary to initialize dataframe
    merged_dict = {}
    merged_dict["bcr_patient_barcode"] = []

    report_path_ex = list(sample_dict.values())[0]
    pathways = list(pd.read_csv(report_path_ex)["Term"])
    for p in pathways:
        merged_dict[p] = []
    for ssample in list(sample_dict):
        report_path = sample_dict[ssample]
        merged_dict["bcr_patient_barcode"].append(ssample)
        report_df = pd.read_csv(report_path)
        for idx, row in report_df.iterrows():
            cur_pathway = row["Term"]
            # normalized expression of each pathway (gene set)
            nes_val = row["NES"]
            merged_dict[cur_pathway].append(nes_val)

    merged_df = pd.DataFrame().from_dict(merged_dict)
    merged_df.to_csv(ssgea_merged_path)
    return merged_df


# subdir = sys.argv[1]
indir = "transcriptome_raw"
# gene_set = sys.argv[2]
# geneset = "AUNG_GASTRIC_CANCER.v2023.2.Hs.gmt"
geneset = "kegg_hsa.gmt"
# outdir = sys.argv[3]
outdir = "ssgsea_output"

# ssgsea_analysis(geneset, indir, outdir)

# create merged file where index is bcr_patient_barcode
ssgea_merged_fname = "ssgea_merged.csv"
# merges and saves to working directory
ssgea_merged = ssgsea_merge(ssgea_merged_fname, outdir)

"""
Rscript to merge with clinical and MATH data (in numerical cols)
combine with J's code. Write gsea_processing.py in R?
merged_df <- merge(clinical, merged_sample_summary, by = "bcr_patient_barcode")
"""
