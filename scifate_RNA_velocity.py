import scvelo as scv
scv.logging.print_version()
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
    
df_cell_file = sys.argv[1]
full_RNA_folder = sys.argv[2]
new_RNA_folder = sys.argv[3]
output_file = sys.argv[4]

df_cell = pd.read_csv(df_cell_file)
#scv.settings.set_figure_params('scvelo')
input_folder = full_RNA_folder
adata = sc.read(input_folder + 'gene_count.mtx').transpose()
df_cell = pd.read_csv(input_folder + 'df_cell.tsv', delimiter="\t")
df_gene = pd.read_csv(input_folder + 'df_gene.tsv', delimiter="\t")
df_cell.index = df_cell["sample"]
df_gene.index = df_gene["gene_id"]
adata.obs = df_cell
adata.var = df_gene

adata_all = adata
input_folder = new_RNA_folder
adata = sc.read(input_folder + 'gene_count.mtx').transpose()
df_cell = pd.read_csv(input_folder + 'df_cell.tsv', delimiter="\t")
df_gene = pd.read_csv(input_folder + 'df_gene.tsv', delimiter="\t")
df_cell.index = df_cell["sample"]
df_gene.index = df_gene["gene_id"]
adata.obs = df_cell
adata.var = df_gene
adata_new = adata

adata_all.layers['spliced'] = adata_all.X
adata_all.layers['unspliced'] = adata_new.X
adata = adata_all

df_cell = pd.read_csv(df_cell_file)

adata_filter = adata[df_cell["sample"], ]
df_cell.index = df_cell["sample"]
adata_filter.obs = df_cell
adata_ori = adata

adata = adata_filter

scv.pp.filter_and_normalize(adata, min_counts=20, min_counts_u=10, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

transition_probs = scv.tl.transition_matrix(adata)

scipy.io.mmwrite(target = output_file, a = transition_probs)
