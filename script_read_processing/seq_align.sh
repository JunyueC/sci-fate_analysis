input_folder=$1
sample=$2
output_folder=$3
sam_tsv="/net/shendure/vol1/home/cao1025/Download/jvarkit/dist/sam2tsv.jar"
fa_file="/net/shendure/vol1/home/cao1025/reference/fasta/hs37d5_mm10.fa"

echo generate alignment file: $sample
java -jar $sam_tsv -r $fa_file $input_folder/$sample.sam >$output_folder/$sample.align