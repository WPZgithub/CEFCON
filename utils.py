import numpy as np
import pandas as pd
from scipy import sparse
import networkx as nx
import scanpy as sc
from typing import Optional, Union
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
import seaborn as sns


def data_preparation(input_expData: Union[str, sc.AnnData, pd.DataFrame],
                     input_priorNet: Union[str, pd.DataFrame],
                     genes_DE: Optional[Union[str, pd.DataFrame, pd.Series]] = None,
                     TF_list: Optional[Union[str, pd.DataFrame, pd.Series]] = None,
                     additional_edges_pct: float = 0.01,
                     only_directed: bool = False) -> sc.AnnData:
    ## Single-cell RNA-seq data
    if isinstance(input_expData, str):
        p = Path(input_expData)
        if p.suffix == '.csv':
            adata = sc.read_csv(input_expData, first_column_names=True)
        else:  # h5ad
            adata = sc.read_h5ad(input_expData)
    elif isinstance(input_expData, sc.AnnData):
        adata = input_expData
    elif isinstance(input_expData, pd.DataFrame):
        adata = sc.AnnData(X=input_expData)
    else:
        raise Exception("Invalid input! The input format must be '.csv' file or '.h5ad' "
                        "formatted hdf5 file, or an 'AnnData' object!", input_expData)
    adata.var_names = adata.var_names.str.upper()
    scData_genes = adata.var_names.values

    ## Prior network data
    if isinstance(input_priorNet, str):
        netData = pd.read_csv(input_priorNet, index_col=None, header=0)
    elif isinstance(input_priorNet, pd.DataFrame):
        netData = input_priorNet
    else:
        raise Exception("Invalid input!", input_priorNet)
    # make sure the genes of prior network are in the input scRNA-seq data
    netData['from'] = netData['from'].str.upper()
    netData['to'] = netData['to'].str.upper()
    netData = netData.loc[netData['from'].isin(scData_genes)
                          & netData['to'].isin(scData_genes), :]
    if only_directed:
        netData = netData.loc[netData['edge_type'] == 'directed', ['from', 'to']]
    else:
        netData = netData.loc[:, ['from', 'to']].drop_duplicates()

    # Transfer into networkx format
    priori_network = nx.from_pandas_edgelist(netData, source='from', target='to', create_using=nx.DiGraph)
    priori_network_nodes = np.array(priori_network.nodes())
    adata.uns['priori_network'] = priori_network

    # Only keep the genes that exist in both single cell data and the reference network
    adata = adata[:, priori_network_nodes]
    if isinstance(adata.X, sparse.csr_matrix):
        node_feature = pd.DataFrame(adata.X.A.T, index=priori_network_nodes)
    else:
        node_feature = pd.DataFrame(adata.X.T, index=priori_network_nodes)

    # in_degree, out_degree and eigenvector (centrality: [0,1])
    in_degree = pd.DataFrame.from_dict(nx.in_degree_centrality(priori_network),
                                       orient='index', columns=['in_degree'])
    out_degree = pd.DataFrame.from_dict(nx.out_degree_centrality(priori_network),
                                        orient='index', columns=['out_degree'])
    # in_eigenvector = pd.DataFrame.from_dict(nx.eigenvector_centrality(priori_network),
    #                                         orient='index', columns=['in_eigenvector'])
    # out_eigenvector = pd.DataFrame.from_dict(nx.eigenvector_centrality(priori_network.reverse()),
    #                                         orient='index', columns=['out_eigenvector'])
    centrality = pd.concat([in_degree, out_degree], axis=1)
    centrality = np.array(centrality.loc[priori_network_nodes, :])
    adata.varm['centrality'] = centrality

    # A map for node index and gene name
    idx_GeneName_map = pd.DataFrame({'idx': range(len(priori_network_nodes)),
                                     'geneName': priori_network_nodes},
                                    index=priori_network_nodes)
    adata.varm['idx_GeneName_map'] = idx_GeneName_map

    edgelist = pd.DataFrame({'from': idx_GeneName_map.loc[netData['from'].tolist(), 'idx'].tolist(),
                             'to': idx_GeneName_map.loc[netData['to'].tolist(), 'idx'].tolist()})

    ## Additional edges with high spearman correlation
    ori_edgeNum = len(edgelist)
    edges_corr = np.absolute(np.array(node_feature.T.corr('spearman')))
    np.fill_diagonal(edges_corr, 0.0)
    x, y = np.where(edges_corr > 0.6)
    addi_top_edges = pd.DataFrame({'from': x, 'to': y, 'weight': edges_corr[x, y]})
    addi_top_k = int(node_feature.shape[0] * (node_feature.shape[0] - 1) * additional_edges_pct)
    if len(addi_top_edges) > addi_top_k:
        addi_top_edges = addi_top_edges.sort_values(by=['weight'], ascending=False)
        addi_top_edges = addi_top_edges.iloc[0:addi_top_k, 0:2]
    edgelist = pd.concat([edgelist, addi_top_edges.iloc[:, 0:2]], ignore_index=True)
    edgelist = edgelist.drop_duplicates(subset=['from', 'to'], keep='first', inplace=False)
    print('{} extra edges (Spearman correlation > 0.6) are added into the prior gene interaction network.'.format(len(edgelist) - ori_edgeNum))
    adata.uns['edgelist'] = edgelist

    ## Differential expression scores
    if genes_DE is not None:
        if isinstance(genes_DE, str):
            genes_DE = pd.read_csv(genes_DE, index_col=0, header=0)
        genes_DE = pd.DataFrame(genes_DE).iloc[:, 0]
        genes_DE.index = genes_DE.index.str.upper()
        genes_DE = genes_DE[genes_DE.index.isin(priori_network_nodes)].abs()
        node_score_auxiliary = pd.Series(np.zeros(len(priori_network_nodes)), index=priori_network_nodes)
        node_score_auxiliary[genes_DE.index] = genes_DE.values
        node_score_auxiliary = np.array(node_score_auxiliary)
        adata.var['node_score_auxiliary'] = node_score_auxiliary

    ## TF list
    is_TF = np.ones(len(priori_network_nodes), dtype=int)
    if TF_list is not None:
        if isinstance(TF_list, str):
            TF_list = pd.read_csv(TF_list, header=None)
        TF_list = pd.DataFrame(TF_list).iloc[:, 0].str.upper()
        is_TF[~np.isin(priori_network_nodes, TF_list)] = 0
    adata.var['is_TF'] = is_TF

    return adata


def gene_clusters(gene_embeddings: pd.DataFrame, resolution: float=1, vis: bool=False, info=None):
    '''
    Clustering genes with provided embeddings using the Leiden community detection
    algorithm in low-dimensional UMAP space.
    '''
    # Clustering genes with obtained embeddings
    adata_gene = sc.AnnData(X=gene_embeddings)
    sc.pp.neighbors(adata_gene, n_neighbors=30, use_rep='X')
    # Higher resolutions lead to more communities, while lower resolutions lead to fewer communities.
    sc.tl.leiden(adata_gene, resolution=resolution) #
    sc.tl.umap(adata_gene, n_components=2, min_dist=0.3)
    pos = adata_gene.obsm['X_umap']
    gene_labels = adata_gene.obs['leiden']

    # draw the graph
    if vis:
        sc.pl.umap(adata_gene, color=['leiden'], legend_loc='on data',
                   legend_fontsize=10, legend_fontoutline=3, title=info)
    return gene_labels, pos, adata_gene


def regulon_activity(ex_matrix: pd.DataFrame, G: nx.DiGraph,
                     out_hub_nodes: Union[set, list], in_hub_nodes:  Union[set, list],
                     DEgenes=None, cell_label=None, num_workers=8):
    '''
    Select top regulons and calculate activity in each cell.
    '''
    from pyscenic.aucell import aucell
    from ctxcore.genesig import GeneSignature

    auc_mtx_out, auc_mtx_in, auc_mtx = None, None, None
    if DEgenes is None:
        DEgenes = set(G.nodes())
    else:
        DEgenes = set(DEgenes)

    # Out_regulons
    out_regulons = {i+'_out({})'.format(len(list(set(G.successors(i)).intersection(DEgenes)))):
                        list(set(G.successors(i)).intersection(DEgenes)) for i in out_hub_nodes
                    if len(set(G.successors(i)).intersection(DEgenes))>=10}
    # In_regulons
    in_regulons = {i+'_in({})'.format(len(list(set(G.predecessors(i)).intersection(DEgenes)))):
                       list(set(G.predecessors(i)).intersection(DEgenes)) for i in in_hub_nodes
                   if len(set(G.predecessors(i)).intersection(DEgenes))>=10}

    # n_cells x n_genes
    out_regulons = [GeneSignature(name=k, gene2weight=v) for k,v in out_regulons.items()]
    in_regulons = [GeneSignature(name=k, gene2weight=v) for k, v in in_regulons.items()]
    regulons = out_regulons + in_regulons
    if len(out_regulons)>0:
        auc_mtx_out = aucell(ex_matrix, out_regulons, num_workers=num_workers, auc_threshold=0.25)
        # Generate a Z-score for each regulon to enable comparison between regulons
        auc_mtx_out_Z = pd.DataFrame(index=auc_mtx_out.index)
        for col in list(auc_mtx_out.columns):
            auc_mtx_out_Z[col] = (auc_mtx_out[col] - auc_mtx_out[col].mean()) / auc_mtx_out[col].std(ddof=0)

    if len(in_regulons) > 0:
        auc_mtx_in = aucell(ex_matrix, in_regulons, num_workers=num_workers, auc_threshold=0.25, normalize=False)
        # Generate a Z-score for each regulon to enable comparison between regulons
        auc_mtx_in_Z = pd.DataFrame(index=auc_mtx_in.index)
        for col in list(auc_mtx_in.columns):
            auc_mtx_in_Z[col] = (auc_mtx_in[col] - auc_mtx_in[col].mean()) / auc_mtx_in[col].std(ddof=0)

    if (len(out_regulons) > 0) & (len(in_regulons) > 0):
        auc_mtx = pd.merge(auc_mtx_out, auc_mtx_in, how='inner', left_index=True, right_index=True)
        # Generate a Z-score for each regulon to enable comparison between regulons
        auc_mtx_Z = pd.DataFrame(index=auc_mtx.index)
        for col in list(auc_mtx.columns):
           auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

    # figure
    # Create a categorical palette to identify the networks
    if cell_label is not None:
        if isinstance(cell_label,pd.DataFrame):
            if cell_label.shape[1]==2:
                network_lable = cell_label.iloc[:, 0]
                network_pal = sns.husl_palette(len(network_lable.unique()), h=.5)
                network_lut = dict(zip(map(str, network_lable.unique()), network_pal))
                network_color = pd.Series(list(network_lable), index=auc_mtx.index, name='CellTypes').map(network_lut)

                network_lable_2 = cell_label.iloc[:, 1]
                network_pal_2 = sns.color_palette('Accent', len(network_lable_2.unique()))
                network_lut_2 = dict(zip(map(str, network_lable_2.unique()), network_pal_2))
                network_color_2 = pd.Series(list(network_lable_2), index=auc_mtx.index, name='CellClusters').map(network_lut_2)
                network_colors = pd.DataFrame(network_color).join(pd.DataFrame(network_color_2))
        else:
            network_lable = cell_label
            network_pal = sns.husl_palette(len(network_lable.unique()), h=.5)
            network_lut = dict(zip(map(str, network_lable.unique()), network_pal))
            network_colors = pd.Series(list(network_lable), index=auc_mtx.index).map(network_lut)

        # in+out
        plt.figure(figsize=(9, 21), dpi=300)
        sns.set_theme(font_scale=1.5)
        g = sns.clustermap(auc_mtx.T, method='ward', square=False, linecolor='black', z_score=0, vmin=-2.5, vmax=2.5,
                           col_cluster=False, col_colors=network_colors, cmap="YlGnBu", figsize=(9, 21),
                           xticklabels=False, yticklabels=True) #cmap="YlGnBu"
        g.cax.set_visible(True)
        g.ax_heatmap.set_ylabel('Regulon-like gene modules')
        g.ax_heatmap.set_xlabel('Cells')
        for label in network_lable.unique():
            g.ax_col_dendrogram.bar(0, 0, color=network_lut[label], label=label, linewidth=0)
        g.ax_col_dendrogram.legend(title='Cell types', loc="upper left", ncol=1,
                                   bbox_to_anchor=(1.04, 0.90), facecolor='white')
        if isinstance(cell_label,pd.DataFrame):
            if cell_label.shape[1]==2:
                xx = []
                for label in network_lable_2.unique():
                    x = g.ax_col_dendrogram.bar(0, 0, color=network_lut_2[label], label=label, linewidth=0)
                    xx.append(x)
                legend2 = plt.legend(xx, network_lable_2.unique(), loc="upper left", title='Cell clusters',
                                     bbox_to_anchor=(0.48, 0.97), bbox_transform=gcf().transFigure, facecolor='white')
                plt.gca().add_artist(legend2)
        g.ax_row_dendrogram.set_visible(False)
        #plt.title('AUCell clustermap for all RGMs', fontdict={'fontsize': 30}, loc='right')
        plt.show()

        #out / in
        plt.figure(figsize=(9, 16), dpi=300)
        sns.set_theme(font_scale=1.5)
        g1 = sns.clustermap(auc_mtx_out.T, method='ward', square=False, linecolor='black', z_score=0, vmin=-2.5, vmax=2.5,
                            row_cluster=True, col_colors=network_colors, cmap="RdBu_r", figsize=(9, 16),
                            xticklabels=False, yticklabels=True, col_cluster=False) #cmap="YlGnBu"
        g1.cax.set_visible(True)
        g1.ax_heatmap.set_ylabel('Regulon-like gene modules')
        g1.ax_heatmap.set_xlabel('Cells')
        for label in network_lable.unique():
            g1.ax_col_dendrogram.bar(0, 0, color=network_lut[label], label=label, linewidth=0)
        g1.ax_col_dendrogram.legend(title='Cell types', loc="upper left", ncol=1,
                                    bbox_to_anchor=(1.10, 0.90), facecolor='white')
        if isinstance(cell_label,pd.DataFrame):
            if cell_label.shape[1] == 2:
                xx = []
                for label in network_lable_2.unique():
                    x = g.ax_col_dendrogram.bar(0, 0, color=network_lut_2[label], label=label, linewidth=0)
                    xx.append(x)
                legend2 = plt.legend(xx, network_lable_2.unique(), loc="upper left", title='Cell clusters',
                                     bbox_to_anchor=(0.48, 0.97), bbox_transform=gcf().transFigure, facecolor='white')
                plt.gca().add_artist(legend2)
        g1.ax_row_dendrogram.set_visible(False)
        #g1.ax_heatmap.set_title('AUCell clustermap for outgoing RGMs', fontdict={'fontsize': 15})
        plt.show()

        #
        plt.figure(figsize=(9, 16), dpi=300)
        sns.set_theme(font_scale=1.5)
        g2 = sns.clustermap(auc_mtx_in.T, method='ward', square=False, linecolor='black', z_score=0, vmin=-2.5, vmax=2.5,
                            row_cluster=True, col_colors=network_colors, cmap="RdBu_r", figsize=(9, 16),
                            xticklabels=False, yticklabels=True, col_cluster=False) #cmap="YlGnBu"
        g2.cax.set_visible(True)
        g2.ax_heatmap.set_ylabel('Regulon-like gene modules')
        g2.ax_heatmap.set_xlabel('Cells')
        for label in network_lable.unique():
            g2.ax_col_dendrogram.bar(0, 0, color=network_lut[label], label=label, linewidth=0)
        g2.ax_col_dendrogram.legend(title='Cell types', loc="upper left", ncol=1,
                                    bbox_to_anchor=(1.10, 0.90), facecolor='white')
        if isinstance(cell_label,pd.DataFrame):
            if cell_label.shape[1] == 2:
                xx = []
                for label in network_lable_2.unique():
                    x = g.ax_col_dendrogram.bar(0, 0, color=network_lut_2[label], label=label, linewidth=0)
                    xx.append(x)
                legend2 = plt.legend(xx, network_lable_2.unique(), loc="upper left", title='Cell clusters',
                                     bbox_to_anchor=(0.48, 0.97), bbox_transform=gcf().transFigure, facecolor='white')
                plt.gca().add_artist(legend2)
        g2.ax_row_dendrogram.set_visible(False)
        #g1.ax_heatmap.set_title('AUCell clustermap for incoming RGMs', fontdict={'fontsize': 15})
        plt.show()

    results_dict = {'regulons':regulons, 'aucell':auc_mtx,
                    'aucell_out':auc_mtx_out, 'aucell_in':auc_mtx_in}
    return results_dict


def regulon_cluster_cell(auc_mtx, true_cell_label, method='ward', k=None):
    '''
    Cluster cells based on regulon activity matrix by using hierarchical clustering.
    '''
    assert method in {'ward', 'complete', 'average', 'single'}
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score, normalized_mutual_info_score, adjusted_rand_score

    auc_mtx_Z = pd.DataFrame(index=auc_mtx.index, columns=list(auc_mtx.columns))
    for row in list(auc_mtx.index):
        auc_mtx_Z.loc[row, :] = (auc_mtx.loc[row, :] - auc_mtx.loc[row, :].mean()) / auc_mtx.loc[row, :].std(ddof=0)

    if k is not None:
        ac = AgglomerativeClustering(n_clusters=k, affinity='euclidean', linkage=method).fit(auc_mtx_Z)
        predicted_cell_label = ac.labels_
        NMIs = normalized_mutual_info_score(true_cell_label, predicted_cell_label)
        ARIs = adjusted_rand_score(true_cell_label, predicted_cell_label)
        Silhouettes = silhouette_score(auc_mtx, predicted_cell_label, metric='euclidean')
        N_clus = k
    else:
        NMIs, ARIs, N_clus, Silhouettes = [], [], [], []
        max_cluster_num = len(set(true_cell_label)) * 2
        for i in range(2, max_cluster_num):
            ac = AgglomerativeClustering(n_clusters=i, affinity='euclidean', linkage=method).fit(auc_mtx_Z)
            out = ac.labels_
            NMIs.append(normalized_mutual_info_score(true_cell_label, out))
            ARIs.append(adjusted_rand_score(true_cell_label, out))
            N_clus.append(i)
            Silhouettes.append(silhouette_score(auc_mtx, out, metric='euclidean'))

    return {'NMIs':NMIs, 'ARIs':ARIs, 'Silhouettes':Silhouettes,'num_clusters':N_clus}
