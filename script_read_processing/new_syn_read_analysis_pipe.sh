
all_output_folder="/net/shendure/vol1/home/cao1025/Projects/nobackup/181022_time_test/output_STAR_para_1/"
control_cell="/net/shendure/vol1/home/cao1025/Projects/processed_data/181023_scitime_test/data/RData/control_cell_id.csv"
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/hg19_mm10/rmchr.gencode.v19.chr_patch_hapl_scaff.annotation.gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"
core=10
script_folder="/net/shendure/vol1/home/cao1025/analysis_script/scRNA_seq_pipe"
reference_fa="/net/shendure/vol1/home/cao1025/reference/fasta/hs37d5_mm10.fa"
varscan="/net/shendure/vol1/home/cao1025/Download/VarScan.v2.3.9.jar"

script_folder="/net/shendure/vol1/home/cao1025/analysis_script/scitimelapse/"
R_script="/net/shendure/vol1/home/cao1025/analysis_script/sci3/sci3_bash_input_ID_output_core.R"
Rscript="/net/shendure/vol1/home/cao1025/bin/Rscript.3.5"
gene_count_script="/net/shendure/vol1/home/cao1025/analysis_script/scRNA_seq_pipe/sciRNAseq_count.py"

#define the bin of python
python_path="/net/shendure/vol1/home/cao1025/anaconda2/bin/"


# Call SNPs from the control files
### Generate a script with: (1) input folder (2) sample list (3) output folder, merge the single cell sam files in the
# input folder, and then call SNPs from the merged sequences
# first generate a tmp folder for all sorted singe cell bam files
sample_list=$control_cell
input_folder=$all_output_folder/sam_splitted/
output_folder=$all_output_folder/combind_control_bam/

mkdir $output_folder
temp_output=$output_folder/tmp_bam
mkdir $temp_output
for sample in $(cat $sample_list); do echo combining $sample; samtools view -bh $input_folder/$sample.sam |samtools sort - -o $temp_output/$sample.bam; done

echo Merging bam files
samtools merge $output_folder/merged.bam $temp_output/*.bam
echo Calling SNPs
samtools mpileup -B -f $reference_fa $output_folder/merged.bam >$output_folder/output.mpileup
java -jar $varscan mpileup2snp $output_folder/output.mpileup --strand-filter 0 >$output_folder/SNP.vcf

echo analysis done; removing temp files....
rm -rf $temp_output
echo All analysis done.

# Generate align file for each single cell sam file
input_folder=$all_output_folder/sam_splitted/
sample_list=$all_output_folder/barcode_samples.txt
output_folder=$all_output_folder/mutation_count
#SNP_VCF="/net/shendure/vol1/home/cao1025/Projects/nobackup/181022_time_test/output/combind_control_bam/SNP.vcf"
mkdir $output_folder
# generate al alignment files for the output
bash_script=$script_folder/seq_align.sh
$Rscript $R_script $bash_script $input_folder $sample_list $output_folder $core
echo analysis done.

# in this script, I will generate a file and count T -> C mutation in each single cell sam file
input_folder=$all_output_folder/mutation_count/
sample_ID=$all_output_folder/barcode_samples.txt
output_folder=$all_output_folder/new_synthesised_reads/
SNP_VCF=$all_output_folder/combind_control_bam/SNP.vcf
# filter the newly synthesised reads for each reads
$Rscript $script_folder/select_newly_synthesised_read.R $input_folder $sample_ID $output_folder $core $SNP_VCF

# in this script, I will generate a file and identify T -> C mutation in each single cell sam file
input_folder=$all_output_folder
sample_list=$all_output_folder/new_synthesised_reads/sample_id.txt
output_folder=$all_output_folder/new_reads_sam/

mkdir $output_folder
# generate al alignment files for the output
bash_script=$script_folder/extract_new_reads.sh
$Rscript $R_script $bash_script $input_folder $sample_list $output_folder $core
echo analysis done.

# Generate gene count matrix for all files
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