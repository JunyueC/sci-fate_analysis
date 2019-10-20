import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import time
import sys
input_dir = sys.argv[1] + "/"

# The raw counts matrix (E) should be a scipy sparse CSC matrix
# with cells as rows and genes as columns

if os.path.isfile(input_dir + '/gene_count.npz'):
    E = scipy.sparse.load_npz(input_dir + '/gene_count.npz')
else:
    E = scipy.io.mmread(input_dir + '/gene_count.mtx').T.tocsc()
    scipy.sparse.save_npz(input_dir + '/gene_count.npz', E, compressed=True)


genes = np.array(scr.load_genes(input_dir + 'df_gene.tsv', delimiter='\t', column=1))

print('Expression matrix shape: {} rows, {} columns'.format(E.shape[0], E.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(E, expected_doublet_rate=0.05)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
scrub.call_doublets(threshold=0.22)
scrub.plot_histogram()
plt.savefig(input_dir + "/hist1.png")
print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
plt.savefig(input_dir + "/umap.png")

score_threshold_logNorm = 0.22
np.sum(scrub.doublet_scores_obs_ > 0.22)
# output the doublet score to another file for analysis with R
np.savetxt(input_dir + "/doublet_score.csv", scrub.doublet_scores_obs_, delimiter=",")
np.savetxt(input_dir + "/stimulated_doublet_score.csv", scrub.doublet_scores_sim_, delimiter=",")