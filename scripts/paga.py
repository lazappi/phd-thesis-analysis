#!/usr/bin/env python3

import os
import json
import numpy
import pandas as pd
import loompy
import anndata
import scanpy as sc

# Show errors (0), warnings (1), info (2) and hints (3)
sc.settings.verbosity = 3

print('Reading clustering parameters...')
params_path = 'output/03-clustering/parameters.json'

with open(params_path) as params_file:
    params = json.load(params_file)

for param in params:
    if param['Parameter'][0] == 'n_pcs':
        n_pcs = param['Value'][0]

    if param['Parameter'][0] == 'knn':
        knn = param['Value'][0]

print('Reading Loom file...')
loom_path = 'data/processed/03-clustered-sel.loom'

col_names = ['Cell', 'Barcode', 'Dataset', 'Sample', 'Cluster']
obs = dict()

with loompy.connect(loom_path) as loom_con:
    X = loom_con.layers[''][()].T

    for col in col_names:
        obs[col] = loom_con.col_attrs[col]

print('Converting to AnnData...')
adata = anndata.AnnData(X=X, obs=obs)

print('Calculating PCs...')
sc.tl.pca(adata, svd_solver='arpack')

print('Calculating neighbour graph...')
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15)

print('Perfoming PAGA...')
sc.tl.paga(adata, groups='Cluster')

print('Calculating cluster graph layout...')
sc.pl.paga(adata, plot=False)

print('Calculating cell graph layout...')
sc.tl.draw_graph(adata, init_pos='paga')

out_dir = 'output/05-paga'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def con2edges(con, names=None, sparse=True):
    print('Converting connectivity matrix to edges...')
    n = con.shape[0]
    edges = pd.DataFrame(columns=['From', 'To', 'Connectivity'])

    for i in range(n):
        for j in range(i + 1, n):
            if names is not None:
                fr = names[i]
                to = names[j]
            else:
                fr = str(i)
                to = str(j)

            connectivity = con[i, j]
            if sparse and connectivity == 0:
                continue

            entry = {'From' : fr, 'To' : to,
                     'Connectivity' : con[i, j]}
            edges = edges.append(entry, ignore_index=True)

    return edges

print('Outputting cluster edges...')
clust_con = adata.uns['paga']['connectivities'].toarray()
clust_edges = con2edges(clust_con)
clust_edges.to_csv(os.path.join(out_dir, 'cluster_edges.csv'),
                   index=False)

print('Outputting cluster tree edges...')
clust_tree_con = adata.uns['paga']['connectivities_tree'].toarray()
clust_tree_edges = con2edges(clust_tree_con)
clust_tree_edges.to_csv(os.path.join(out_dir, 'cluster_tree_edges.csv'),
                        index=False)

print('Outputting cluster embedding...')
clust_embedding = pd.DataFrame(adata.uns['paga']['pos'], columns=['X', 'Y'])
clust_embedding['Cluster'] = range(clust_embedding.shape[0])
clust_embedding = clust_embedding[['Cluster', 'X', 'Y']]

clust_embedding.to_csv(os.path.join(out_dir, 'cluster_embedding.csv'),
                       index=False)

print('Outputting cell edges...')
cells = adata.obs['Cell']
cell_con = adata.uns['neighbors']['connectivities']
cell_edges = pd.DataFrame(columns=['From', 'To', 'Connectivity'])
n_rows = len(cell_con.indptr)

for i in range(len(cell_con.indptr) - 1):
    row_ind = cell_con.indices[cell_con.indptr[i]:cell_con.indptr[i + 1]]
    print(f'\r\tRow {i} of {n_rows}', end='')
    for k, j in enumerate(row_ind):
        if j > i:
            con = cell_con.data[cell_con.indptr[i] + k]
            fr = cells[i]
            to = cells[j]
            entry = {'From' : fr, 'To' : to, 'Connectivity' : con}
            cell_edges = cell_edges.append(entry, ignore_index=True)
print('\n')

cell_edges.to_csv(os.path.join(out_dir, 'cell_edges.csv'),
                  index=False)

print('Outputting cell embedding...')
x = adata.obsm['X_draw_graph_fa'][:, 0]
y = adata.obsm['X_draw_graph_fa'][:, 1]
cell_embedding = pd.DataFrame({'Cell' : cells, 'X' : x, 'Y' : y})
cell_embedding.to_csv(os.path.join(out_dir, 'cell_embedding.csv'),
                      index=False)

print('Done!')

