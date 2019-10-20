# This script is the main script for screening TF ang gene regulatory links based on the covariance between TF expression
# and gene newly synthesis rate (or expression) as well as the enrichment of TF binding-motif/chip-seq peaks near gene's
# promoter
suppressMessages(require(cowplot))
suppressMessages(require(tidyverse))
suppressMessages(require(monocle))
require(stringr)
require(glmnet)
require(BiocParallel)

require(RcisTarget)

# load the helper functions 
load("./linkage_function.RData")

# get TF list from RcisTarget package: 
data(motifAnnotations_hgnc)
TF_list = unique(motifAnnotations_hgnc$TF)

# define the location of the motif reference downloaded from RcisTarget: https://resources.aertslab.org/cistarget/
# for human motif matrix, it can be downloaded from: https://shendure-web.gs.washington.edu/content/members/cao1025/public/nobackup/sci_fate/data/hg19-tss-centered-10kb-7species.mc9nr.feather
motif_ref = "./data/hg19-tss-centered-10kb-7species.mc9nr.feather"
# define the location of chip-seq peak matrix downloaded from https://amp.pharm.mssm.edu/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets
df_gene_TF_link_ENCODE = "./data/df_TF_gene_ENCODE.RData"

# define core number for linkage analysis (core_n_lasso)
core_n_lasso = 5
# define core number for Rcistarget analysis (core_n_filtering)
core_n_filtering = 1

# define output folder
output_folder = "./output/TF_gene_analysis"
dir.create(output_folder)

# load the monocle2 cds object for full gene expression and newly synthesized reads of each cell
# downloaded from: https://shendure-web.gs.washington.edu/content/members/cao1025/public/nobackup/sci_fate/data/cds_all_new.RDS
load("./data/cds_all_new.RDS")

# check if the cell names and gene names are the same
all(rownames(cds_all) == rownames(cds_new))
all(colnames(cds_all) == colnames(cds_new))

# screening TF ang gene regulatory links based on the covariance between TF expression
# and gene newly synthesis rate as well as the enrichment of TF binding-motif/chip-seq peaks near gene's
# promoter

# Main parameters:
# cds_all: monocle2 cds object for full gene expression of cells
# cds_new: monocle2 cds object for newly synthesized gene expression of cells
# TF_list: gene names of TFs for linkage analysis
# output_links_new_RNA: output folder
# gene_filter_rate: minimum percentage of expressed cells for gene filtering 
# cell_filter_UMI: minimum number of UMIs for cell filtering
# core_n_lasso: number of cores for lasso regression in linkage analysis
# core_n_filtering: number of cores for filtering TF-gene links
# motif_ref: TF binding motif data as described above
# df_gene_TF_link_ENCODE: TF chip-seq data as described above

output_links_new_RNA = file.path(output_folder, "Link_new")
dir.create(output_links_new_RNA)
cds_TF_gene_linkage_analysis(cds_all, cds_new, TF_list, output_links_new_RNA,
                                        gene_filter_rate = 0.1, cell_filter_UMI = 10000,
                                        core_n_lasso, core_n_filtering,
                                        motif_ref, df_gene_TF_link_ENCODE)

# output file:
# "./output/TF_gene_analysis/Link_new/df_gene_TF_all.RDS" is the output file and include the following columns:
# TF: target TF
# linked_gene: linked gene for the target TF
# Conf: link validated by binding motif or chip-seq peaks
# linked_gene_id: linked gene id for the target TF
# linked_gene_type: linked gene type for the target TF
# corcoef: correlation coefficient recovered from lasso regression
# r_2: computed r squared for predicting linked gene expression in lasso regression
# TF_link: identified TF and gene links
