#!/bin/bash
# this is the file for generating newly synthesized gene expression matrix for sci-fate

# define the output folder (same as the Main_script_full_expression_gene_count_matrix.sh)
all_output_folder=$1
# define the VCF file for background mutation
SNP_VCF=$2
# define the single cell sample names for generating newly synthesized gene expression matrix
sample_list=$3
sample_ID=$sample_list

# define the gtf file for gene counting
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/hg19_mm10/rmchr.gencode.v19.chr_patch_hapl_scaff.annotation.gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"
# define core number for processing
core=5
# define the genome fasta file for mutation check
reference_fa="/net/shendure/vol1/home/cao1025/reference/fasta/hs37d5_mm10.fa"

# define the location of script folder 
script_folder="./script_read_processing/"
R_script=$script_folder/sci3_bash_input_ID_output_core.R
gene_count_script=$script_folder/sciRNAseq_count.py

#define the bin of python
python_path="/net/shendure/vol1/home/cao1025/anaconda2/bin/"
# define the location of VarScan for mutation count
varscan="/net/shendure/vol1/home/cao1025/Download/VarScan.v2.3.9.jar"
# define the location of R script
Rscript="/net/shendure/vol1/home/cao1025/bin/Rscript.3.5"

# Generate alignment file for each single cell sam file
input_folder=$all_output_folder/sam_splitted/
output_folder=$all_output_folder/mutation_count
mkdir $output_folder
bash_script=$script_folder/seq_align.sh
$Rscript $R_script $bash_script $input_folder $sample_list $output_folder $core
echo analysis done.

# in this script, I will generate a file to count T -> C mutations in each single cell sam file
input_folder=$all_output_folder/mutation_count/
output_folder=$all_output_folder/new_synthesised_reads/
# filter the newly synthesised reads for each reads
$Rscript $script_folder/select_newly_synthesised_read.R $input_folder $sample_ID $output_folder $core $SNP_VCF

# in this script, I will generate a single cell sam file with T -> C mutations for each cell
input_folder=$all_output_folder
sample_list=$all_output_folder/new_synthesised_reads/sample_id.txt
output_folder=$all_output_folder/new_reads_sam/
mkdir $output_folder
# generate al alignment files for the output
bash_script=$script_folder/extract_new_reads.sh
$Rscript $R_script $bash_script $input_folder $sample_list $output_folder $core
echo analysis done.

# Generate gene count matrix for newly synthesized reads for each cell
input_folder=$all_output_folder/new_reads_sam/
sample_ID=$all_output_folder/new_synthesised_reads/sample_id.txt
output_folder=$all_output_folder/report/newly_syn_human_mouse_gene_count/

script=$gene_count_script
echo "Start the gene count...."
$python_path/python $script $gtf_file $input_folder $sample_ID $core

echo "Make the output folder and transfer the files..."
mkdir -p $output_folder
cat $input_folder/*.count > $output_folder/count.MM
rm $input_folder/*.count
cat $input_folder/*.report > $output_folder/report.MM
rm $input_folder/*.report
mv $input_folder/*_annotate.txt $output_folder/
echo "All output files are transferred~"
