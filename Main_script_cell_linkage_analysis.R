# This script is the main script for linking future and past cell states in sci-fate
suppressMessages(require(cowplot))
suppressMessages(require(tidyverse))
suppressMessages(require(monocle))
load("./script_downstream_processing/linkage_function.RData")

# define output folder
output_folder = "./output/cell_linkage_analysis"
dir.create(output_folder)

# load the monocle cds object for full gene expression and newly synthesized reads
# downloaded from: https://shendure-web.gs.washington.edu/content/members/cao1025/public/nobackup/sci_fate/data/cds_all_new.RDS
load("./data/cds_all_new.RDS")

# define cell time labeles of different treatment time
cell_groups = cds_all$treatment_time

# define the control and treatment time label for sci-fate capture rate estimation
control_group = "0h"
stimu_group = "2h"

# define the nearby groups of sci-fate; each group are organized as c(past_time_point, future_time_point)
# the time interval should be equal with the sci-fate labeling time
group_pairs = list(c("0h", "0h"), c("0h", "2h"), 
                  c("2h", "4h"),
                  c("4h", "6h"),
                  c("6h", "8h"),
                  c("8h", "10h"))

# define core number for analysis
core_num = 3

# check if the cell names of full expression matrix and newly synthesized gene matrix are the same
all(rownames(cds_all) == rownames(cds_new))
all(colnames(cds_all) == colnames(cds_new))

# estimate the sci-fate capture rate
cat("\nEstimate capture rate...")
Detection_rate = estimate_detection_rate(cds_all, cds_new, cell_groups, control_group, stimu_group, output_folder)

# estimate the gene degradation rate
cat("\nEstimate gene degradation rate...")
df_gene_degradation = estimate_degradation_rate(cds_all, cds_new, cell_groups, control_group, stimu_group, group_pairs, output_folder, Detection_rate)

# for each cell compute its past transcriptome
cat("\nEstimate cell past state...")
cds_past = estimate_past_cell_state(cds_all, cds_new, cell_groups, control_group, stimu_group, group_pairs, output_folder, Detection_rate, df_gene_degradation)

require(Seurat)

# for each group pair, integrate the data to identify the past and future state of each cell
linked_cell_output_folder = file.path(output_folder, "Cell_linkage_anlysis")
if(!file.exists(linked_cell_output_folder)) {
    dir.create(linked_cell_output_folder)
}
df_result = cds_linkage_analysis_compute_distance_correlation(cds_all, cds_past, cds_new, cell_groups, group_pairs, linked_cell_output_folder, df_gene_degradation, core_num)
output_folder_transition_prob = file.path(output_folder, "Transition_prob")
cds_compute_transition_prob(cds_all, cds_new, cell_groups, group_pairs, output_folder_transition_prob)
df_combined = combine_distance_cor_prob(linked_cell_output_folder, output_folder_transition_prob)
df_top_links = select_top_links(df_combined, cell_groups, group_pairs)
df_top_links = df_top_links %>% filter(group != "0h_0h")
saveRDS(df_top_links, file = file.path(output_folder, "df_cell_linkages.RDS"))

# Construct a single cell trajectory for each cell
df_cell = readRDS("./data/df_main_umap.RDS")
df_cell$sample = str_replace_all(df_cell$sample, pattern = "-", replacement = ".")
df_cell$umap_1 = df_cell$umap_1_combined
df_cell$umap_2 = df_cell$umap_2_combined
df_trajectory = find_trajectory_all_cells(df_cell, df_top_links, core_num)
df_output = df_trajectory %>% select(cell_id, time_point = treatment_time_num, trajectory_cell_name = sample)
saveRDS(df_output, file = file.path(output_folder, "df_trajectory.RDS"))

# output: "./output/cell_linkage_analysis/df_trajectory.RDS" is a R RDS file
# cell_id: cell id of each trajectory
# time_point: time point in the single cell trajectory
# trajectory_cell_name: linked cell in the target time point within the single cell trajectory
