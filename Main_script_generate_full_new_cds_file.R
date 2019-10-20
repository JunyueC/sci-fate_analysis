# this script accept the all input folder, and then generate cds file for full expression matrix and newly synthesized
# gene expression matrix

suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(require(monocle))

# load the helper functions 
load("./script_downstream_processing/linkage_function.RData")

# all output folder used in Main_script_full_expression_gene_count_matrix.sh or Main_script_new_expression_gene_count_matrix.sh
all_output_folder = "./all_output_folder"


# construct monocle2 cds file for the full gene expression
report_folder = file.path(all_output_folder, "report", "human_mouse_gene_count")
result = sciRNAseq_gene_count_summary(report_folder)
cds_all = cds_construct(result[[3]], result[[1]], result[[2]])

# construct monocle2 cds file for the newly synthesized gene expression
report_folder = file.path(all_output_folder, "report", "newly_syn_human_mouse_gene_count")
result = sciRNAseq_gene_count_summary(report_folder)
cds_new = cds_construct(result[[3]], result[[1]], result[[2]])

# Save the cds_all and cds_new into the output folder
save(cds_all, cds_new, file = file.path(all_output_folder, "cds_all_new.RData"))
