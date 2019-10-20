input_folder=$1 # input folder should include both sam_splitted folder and new_synthesised_reads folder
sample=$2
output_folder=$3

input_sam=$input_folder/sam_splitted/$sample.sam
input_read=$input_folder/new_synthesised_reads/$sample.txt
picard_source="/net/shendure/vol1/home/cao1025/Download/picard.jar"
mkdir $output_folder
java -jar $picard_source FilterSamReads I=$input_sam O=$output_folder/$sample.sam READ_LIST_FILE=$input_read FILTER=includeReadList
echo analyis done: $sample