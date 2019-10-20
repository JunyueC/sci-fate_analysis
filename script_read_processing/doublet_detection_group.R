#!/net/shendure/vol1/home/cao1025/bin/Rscript.3.5
suppressMessages(library(monocle))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
load("~/Projects/notebook/sciRNAseq/scRNA.image.RData")

cds_file = "/net/shendure/vol1/home/cao1025/Projects/processed_data/181103_sciMM/data/Main_A549/181211_cds_A549_filtered.RDS"
doublet_output_folder="/net/shendure/vol1/home/cao1025/Projects/processed_data/181103_sciMM/data/Main_A549/doublet_analysis_group/"

cds = readRDS(cds_file)
groups = cds$RT_group

# Here I define a function, that accept a cds, a output folder, and then perform doublet analysis with the cds,
# and return a data frame including the doublet score in the output folder
cds_doublet_analysis <- function(cds, output_folder, python_home="/net/shendure/vol1/home/cao1025/anaconda3/bin/python3", script="/net/shendure/vol1/home/cao1025/analysis_script/scitimelapse/doublet_detector.py") {
    cat("\nOutput cds into the output folder....")
    cds_to_matrix(cds, output_folder)
    
    cat("\n*******Start doublet analysis*******")
    cat("\npython home: ", python_home)
    cat("\npython script: ", script)
    system(paste(python_home, script, output_folder, sep = " "), intern = T, wait = TRUE)
    cat("\n Doublet analysis complete")
    
    cat("\nSummarizing the data...")
    input_folder=output_folder
    cells = read_csv(paste(input_folder, "df_cell.tsv", sep = "/"))
    doublet_score = read_csv(paste(input_folder, "doublet_score.csv", sep = "/"), col_names = F)
    stimu_doublet_score = read_csv(paste(input_folder, "stimulated_doublet_score.csv", sep = "/"), col_names = F)
    cells$db_score = doublet_score$X1
    df_cell_db = cells
    df_stim_db = stimu_doublet_score
    
    #g1 = (ggplot()
    #  + geom_histogram(aes(x=df_stim_db$X1, fill = "stimulated data"), alpha = 0.5)
    #  + geom_histogram(aes(x=df_cell_db$db_score, fill = "Real data"), alpha = 0.5)
    #  + labs(x="Doublet score", y = "Cell count")
    #  + scale_fill_discrete(""))
    #write_csv(df_cell_db, file.path(input_folder, "df_cell_dbscore.csv"))
    #write_csv(df_stim_db, file.path(input_folder, "df_cell_stimulated_doublet_score.csv"))
    #cowplot::save_plot(g1, filename = paste(input_folder, "hist_doublet.pdf", sep = "/"))
    
    cat("\nAll analysis done")
    return(list(df_cell_db, df_stim_db))
}


df_cell = NULL
df_stimu = NULL

# perfrom doublet analysis for each group
unique_group = unique(groups)
for(x in unique_group) {
    cat("\nProcessing group: ", x)
    tmp_output_folder = file.path(doublet_output_folder, as.character(x))
    tmp_output = cds_doublet_analysis(cds[, groups == x], tmp_output_folder)
    tmp_df_cell = tmp_output[[1]]
    tmp_df_stim = tmp_output[[2]]
    tmp_df_cell$doublet_group = x
    tmp_df_stim$doublet_group = x
    df_cell = rbind(df_cell, tmp_df_cell)
    df_stimu = rbind(df_stimu, tmp_df_stim)
    cat("\nCleaning intermediate files: ", x)
    unlink(tmp_output_folder, recursive=TRUE)
}

cat("\noutput the files into the output folder: ", doublet_output_folder)
write_csv(df_cell, file.path(doublet_output_folder, "df_cell_dbscore.csv"))
write_csv(df_stimu, file.path(doublet_output_folder, "df_cell_stimulated_doublet_score.csv"))

g1 = (ggplot()
      + geom_histogram(aes(x=df_stimu$X1, fill = "stimulated data"), alpha = 0.5)
      + geom_histogram(aes(x=df_cell$db_score, fill = "Real data"), alpha = 0.5)
      + labs(x="Doublet score", y = "Cell count")
      + scale_fill_discrete(""))

cowplot::save_plot(g1, filename = paste(doublet_output_folder, "hist_doublet.pdf", sep = "/"))

# Here I define a function that accept a vector of real cell doublet score, and a vector of doublet score for stimulated cells
# and a threshold, and it print out the detected doublet rate in real data and stimulated data and estimated doublet score
# for the whole data set
cal_doublet_rate <- function(real_score, sti_score, thresh = 0.25) {
    read_db = real_score
    sti_db = sti_score
    thresh_value = thresh
    cat("\nPercent of real data detected as doublets: ", sum(read_db > thresh_value) / length(read_db))
    cat("\nPercent of stimulated data detected as doublets: ", sum(sti_db > thresh_value) / length(sti_db))
    cat("\nPredicted doublet ratio in the real data: ", sum(read_db > thresh_value) / length(read_db) / (sum(sti_db > thresh_value) / length(sti_db)))
}


cat("\n*******Using 0.15 as threshold: ")
cal_doublet_rate(df_cell$db_score, df_stimu$X1, 0.15)

cat("\n*******Using 0.18 as threshold: ")
cal_doublet_rate(df_cell$db_score, df_stimu$X1, 0.18)

cat("\n*******Using 0.2 as threshold: ")
cal_doublet_rate(df_cell$db_score, df_stimu$X1, 0.2)

cat("\n*******Using 0.22 as threshold: ")
cal_doublet_rate(df_cell$db_score, df_stimu$X1, 0.22)

cat("\n*******Using 0.25 as threshold: ")
cal_doublet_rate(df_cell$db_score, df_stimu$X1, 0.25)


cat("\nAll analysis done")