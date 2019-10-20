# sci-fate_analysis

This folder include the main script for (1) Generating single cell gene count matrix for both full RNA expression and newly synthesized gene expression (2) link TF and gene targets based on the covariance between TF expression and gene synthesis rate and TF binding data (3) infer single cell state transition trajectory in sci-fate.

The function of main scripts:

Main_script_full_expression_gene_count_matrix.sh:
Accept a input fastq folder of sci-RNA-seq, and then use the default parameter for generating single cell gene count matrix of full RNA expression

Main_script_new_expression_gene_count_matrix.sh:
Generate single cell gene expression matrix of newly synthesized reads for sci-fate

Main_script_generate_full_new_cds_file.R:
Generate monocle cds file for full expression matrix and newly synthesized gene expression matrix for downstream analysis

Main_script_TF_gene_linkage.R:
This script is the main script for screening TF ang gene regulatory links based on the covariance between TF expression and gene newly synthesis rate as well as the enrichment of TF binding-motif/chip-seq peaks near gene's promoter.

Main_script_cell_linkage_analysis.R:
This script is the main script for linking future and past cell states in sci-fate
