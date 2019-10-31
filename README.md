# sci-fate_analysis

This folder includes the main scripts for (1) Generating single cell gene count matrix for both full mRNA expression and newly synthesized gene expression; (2) Link TF and gene targets based on the covariance between TF expression and gene synthesis rate as well as TF binding data; (3) Infer single cell state transition trajectories as shown in the sci-fate manuscript. 

The function of main scripts:

Main_script_full_expression_gene_count_matrix.sh:
Generate single cell gene count matrix of full mRNA expression.

Main_script_new_expression_gene_count_matrix.sh:
Generate single cell gene expression matrix of newly synthesized mRNA expression for sci-fate.

Main_script_generate_full_new_cds_file.R:
Generate Monocle cds file for full gene expression matrix and newly synthesized gene expression matrix used for downstream analysis.

Main_script_TF_gene_linkage.R:
Screen TF and gene regulatory links based on the covariance between TF expression and gene newly synthesis rate as well as the enrichment of TF binding-motif/chip-seq peaks near genes' promoters.

Main_script_cell_linkage_analysis.R:
Construct single cell state transition trajectory in sci-fate. 

Some of the main scripts depends on R package Monocle3/alpha version. Its source code and installation instructions are both included in folder: Monocle_Package_source.
