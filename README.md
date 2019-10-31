# sci-fate_analysis

This folder include the main scripts for (1) Generating single cell gene count matrix for both full RNA expression and newly synthesized gene expression (2) Link TF and gene targets based on the covariance between TF expression and gene synthesis rate as well as TF binding data (3) Infer single cell state transition trajectories as shown in sci-fate manuscript. 

The function of main scripts:

Main_script_full_expression_gene_count_matrix.sh:
Generate single cell gene count matrix of full RNA expression

Main_script_new_expression_gene_count_matrix.sh:
Generate single cell gene expression matrix of newly synthesized reads for sci-fate

Main_script_generate_full_new_cds_file.R:
Generate monocle cds file for full expression matrix and newly synthesized gene expression matrix for downstream analysis

Main_script_TF_gene_linkage.R:
Screen TF and gene regulatory links based on the covariance between TF expression and gene newly synthesis rate as well as the enrichment of TF binding-motif/chip-seq peaks near gene's promoter.

Main_script_cell_linkage_analysis.R:
the main script for constructing single cell state transition trajectory in sci-fate. 

As the analysis is done with monocle3/alpha version, we also include its source code and installation instructions in folder: Monocle_Package_source.
