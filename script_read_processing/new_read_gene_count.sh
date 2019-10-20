
all_output_folder="/net/shendure/vol1/home/cao1025/Projects/nobackup/181022_time_test/output/"
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/hg19_mm10/rmchr.gencode.v19.chr_patch_hapl_scaff.annotation.gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"
core=5
script_folder="/net/shendure/vol1/home/cao1025/analysis_script/scRNA_seq_pipe"

#define the bin of python
python_path="/net/shendure/vol1/home/cao1025/anaconda2/bin/"

input_folder=$all_output_folder/new_reads_sam/
sample_ID=$all_output_folder/barcode_samples.txt

################# gene count
# count reads mapping to genes
output_folder=$all_output_folder/report/newly_syn_human_mouse_gene_count/

script=$script_folder/sciRNAseq_count.py
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