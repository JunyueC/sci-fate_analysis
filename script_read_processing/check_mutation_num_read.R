
require("tidyverse")
require("BiocParallel")

#input_folder="~/Projects/nobackup/181022_time_test/output/mutation_count/"
#sample_ID="~/Projects/nobackup/181022_time_test/output/barcode_samples.txt"
#output_folder = "~/Projects/nobackup/181022_time_test/output/new_synthesised_reads"
#core = 10
#SNP_VCF = "~/Projects/nobackup/181022_time_test/output/combind_control_bam/SNP.vcf"

args = commandArgs(trailingOnly=TRUE)

input_folder=args[1]
sample_ID=args[2]
output_folder=args[3]
core = as.numeric(args[4])
SNP_VCF = args[5]

quality_filter = 45
end_dist = 0
target_mut_filter_rate = 0.3

if(!file.exists(output_folder)) {
    dir.create(output_folder)
}

sample_names = (read_csv(sample_ID, col_names = F))$X1

#sample_names = sample_names[1:5]

SNP = read.table(SNP_VCF, sep =  "\t", header = T)
SNP$chr_pos = str_c(SNP$Chrom, SNP$Position, SNP$Ref, SNP$Var, sep = "-")

function_sample <- function(target_id, input_folder, output_folder, core, SNP, quality_filter = 45, end_dist = 0, target_mut_filter_rate = 0.3) {
    cat("Process sample: ", target_id)
    cat("\n")
    align_file = file.path(input_folder, paste0(target_id, ".align"))
    test_input = read_tsv(align_file, col_types = cols(.default = "c"))
    colnames(test_input)[1] = "READ_NAME"
    
    ori_test = test_input
    ori_test$READ_POS = as.numeric(as.character(ori_test$READ_POS))
    ori_test$REF_POS = as.numeric(as.character(ori_test$REF_POS))
    ori_test$FLAG = as.numeric(as.character(ori_test$FLAG))

    test_input = ori_test
    test_input = test_input %>% filter(!is.na(CHROM))
    test_input = test_input %>% filter(BASE != ".", REF != ".")
    test_input$REF = str_to_upper(test_input$REF)
    test_input$BASE = str_to_upper(test_input$BASE)
    test_input = test_input %>% filter((BASE) != (REF))
    if(nrow(test_input) == 0) {
        return(-1)
    }
    # test_input  = test_input %>% filter(QUAL %in% c("A", "E"))

    #length(unique((test_input %>% filter((FLAG == 0 & REF == "T" & BASE == "C") | (FLAG == 16 & REF == "A" & BASE == "G")))$READ_NAME)) / length(unique(ori_test$READ_NAME)) 
    # Remove SNPs
    test_input$chr_pos = str_c(test_input$CHROM, test_input$REF_POS, test_input$REF, test_input$BASE, sep = "-")
    test_input = test_input %>% filter(!(chr_pos %in% SNP$chr_pos))
    
    if(nrow(test_input) == 0) {
        return(-1)
    }
    test_input = test_input[(sapply(test_input$QUAL, utf8ToInt)) > quality_filter, ]
    if(nrow(test_input) == 0) {
        return(-1)
    }
    #length(unique((test_input %>% filter((FLAG == 0 & REF == "T" & BASE == "C") | (FLAG == 16 & REF == "A" & BASE == "G")))$READ_NAME)) / length(unique(ori_test$READ_NAME)) 
    # remove based on mutation distance to end point
    end_pos = ori_test %>% filter(!is.na(READ_POS)) %>% group_by(READ_NAME) %>% summarise(end_point = max(READ_POS))
    test_input = left_join(test_input, end_pos %>% select(READ_NAME, end_point))
    tmp = test_input %>% filter(READ_POS > end_dist, READ_POS < end_point - end_dist)
    if(nrow(tmp) == 0) {
        return(-1)
    }
    #length(unique((tmp %>% filter((FLAG == 0 & REF == "T" & BASE == "C") | (FLAG == 16 & REF == "A" & BASE == "G")))$READ_NAME)) / length(unique(ori_test$READ_NAME)) 

    # remove distance based on the ratio of T-> C mutation
    tmp_mut_num = tmp %>% group_by(READ_NAME) %>% summarise(mut_num = n())
    tmp_mut_num = tmp %>% group_by(READ_NAME) %>% summarise(mut_num = n())
    tmp_target_mut_num = (tmp %>% filter((FLAG == 0 & REF == "T" & BASE == "C") | (FLAG == 16 & REF == "A" & BASE == "G"))
                          %>% group_by(READ_NAME) %>% summarise(target_mut_num = n()))
    tmp_target_mut_num = tmp_target_mut_num %>% select(target_mut_num)
    write_csv(tmp_target_mut_num, file.path(output_folder, paste0(target_id, ".csv")), col_names = F)
    
    return(0)
}

newly_synthesised_read <- bplapply(sample_names, function(target_id) {
    function_sample(target_id, input_folder = input_folder, output_folder = output_folder, core = core, SNP = SNP, quality_filter = quality_filter, end_dist = end_dist, target_mut_filter_rate = target_mut_filter_rate)
}, BPPARAM = MulticoreParam(workers = core))