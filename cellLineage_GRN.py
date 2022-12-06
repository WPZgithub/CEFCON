from typing import Optional
from torch_geometric.typing import Adj
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import trange
import math
from scanpy import AnnData

import torch
from torch import Tensor, nn
import torch.nn.functional as F

import torch_geometric as pyg
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.nn.inits import glorot, zeros
from torch_geometric.nn import Linear, BatchNorm, DeepGraphInfomax
from torch_geometric.utils import (
    remove_self_loops,
    add_self_loops,
    softmax,
    to_undirected,
)
from torch_geometric.data import Data


class GraphAttention_layer(MessagePassing):
    def __init__(self,
                 input_dim: int,
                 output_dim: int,
                 attention_type: str = 'COS',
                 flow: str = 'source_to_target',
                 heads: int = 1,
                 concat: bool = True,
                 dropout: float = 0.0,
                 add_self_loops: bool = True,
                 to_undirected: bool = False,
                 **kwargs):
        kwargs.setdefault('aggr', 'add')
        kwargs.setdefault('flow', flow)
        super(GraphAttention_layer, self).__init__(node_dim=0, **kwargs)

        self.input_dim = input_dim
        self.output_dim = output_dim
        self.heads = heads
        self.concat = concat
        self.dropout = dropout
        self.add_self_loops = add_self_loops
        self.attention_type = attention_type
        self.to_undirected = to_undirected

        # AD: Additive; SD: abs(scaled-dot product); COS: abs(cosine)
        assert attention_type in ['SD', 'COS', 'AD']
        self.lin_l = Linear(input_dim, heads * output_dim, bias=False,
                            weight_initializer='glorot')
        self.lin_r = Linear(input_dim, heads * output_dim, bias=False,
                            weight_initializer='glorot')

        if self.attention_type == 'AD':
            self.att_l = nn.Parameter(Tensor(1, heads, output_dim))
            self.att_r = nn.Parameter(Tensor(1, heads, output_dim))
        else:
            self.register_parameter('att_l', None)
            self.register_parameter('att_r', None)

        if concat:
            self.bias = nn.Parameter(Tensor(heads * output_dim))
            self.weight_concat = nn.Parameter(Tensor(heads * output_dim, output_dim))
        else:
            self.bias = nn.Parameter(Tensor(output_dim))
            self.register_parameter('weight_concat', None)

        self._alpha = None
        self.reset_parameters()

    def reset_parameters(self):
        self.lin_l.reset_parameters()
        self.lin_r.reset_parameters()
        glorot(self.att_l)
        glorot(self.att_r)
        zeros(self.bias)
        glorot(self.weight_concat)

    def forward(self, x: Tensor, edge_index: Adj, x_auxiliary: Tensor,
                return_attention_weights: Optional[bool] = None):
        N, H, C = x.size(0), self.heads, self.output_dim

        if self.to_undirected:
            edge_index = to_undirected(edge_index)
        if self.add_self_loops:
            edge_index, _ = remove_self_loops(edge_index)
            edge_index, _ = add_self_loops(edge_index, num_nodes=N)
        else:
            edge_index, _ = remove_self_loops(edge_index)

        x_l = self.lin_l(x).view(-1, H, C)
        x_r = self.lin_r(x).view(-1, H, C)

        if self.attention_type == 'AD':
            out = self.propagate(edge_index, x=(x_l, x_r), x_norm=None,
                                 x_auxiliary=x_auxiliary, size=None)
        elif self.attention_type == 'COS':
            x_norm_l = F.normalize(x_l, p=2., dim=-1)
            x_norm_r = F.normalize(x_r, p=2., dim=-1)
            out = self.propagate(edge_index, x=(x_l, x_r), x_norm=(x_norm_l, x_norm_r),
                                 x_auxiliary=x_auxiliary, size=None)
        else:  # SD
            out = self.propagate(edge_index, x=(x_l, x_r), x_norm=None,
                                 x_auxiliary=x_auxiliary, size=None)

        alpha = self._alpha
        self._alpha = None

        if self.concat:
            out = out.view(-1, self.heads * self.output_dim)
            out += self.bias
            out = torch.matmul(out, self.weight_concat)
        else:
            out = out.mean(dim=1)
            out += self.bias

        if isinstance(return_attention_weights, bool):
            assert alpha is not None
            return out, (edge_index, alpha)
        else:
            return out

    def message(self, edge_index_i: Tensor, x_i: Tensor, x_j: Tensor,
                x_norm_i: Optional[Tensor], x_norm_j: Optional[Tensor],
                x_auxiliary_j: Tensor, size_i: Optional[int]):
        T = 1.0
        if self.attention_type == 'AD':
            alpha = (x_j * self.att_l).sum(-1) + (x_i * self.att_r).sum(-1)
            alpha = x_auxiliary_j * F.leaky_relu(alpha, 0.2)
        elif self.attention_type == 'COS':
            alpha = x_auxiliary_j * torch.abs((x_norm_i * x_norm_j).sum(dim=-1))
            T = 0.25 #0.1, 0.25, 0.5, 0.75, 1.0
        else:  # 'SD'
            alpha = x_auxiliary_j * torch.abs((x_i * x_j).sum(dim=-1)) / math.sqrt(self.output_dim)

        alpha = softmax(alpha / T, edge_index_i, num_nodes=size_i) # with Temperature hyperparameter
        self._alpha = alpha
        alpha = F.dropout(alpha, p=self.dropout, training=self.training)
        return x_j * alpha.view(-1, self.heads, 1)

    def __repr__(self):
        return '{}({}, {}, heads={}, type={})'.format(self.__class__.__name__,
                                                      self.input_dim,
                                                      self.output_dim,
                                                      self.heads,
                                                      self.attention_type)


class GRN_Encoder(nn.Module):
    def __init__(self,
                 input_dim: int,
                 hidden_dim: int,
                 output_dim: int,
                 heads_num: int = 1,
                 dropout: float = 0.0,
                 attention_type: str = 'COS'):
        super(GRN_Encoder, self).__init__()
        assert attention_type in ['SD', 'COS', 'AD']

        self.att_weights_first = None
        self.att_weights_second = None
        self.x_embs = None

        self.x_input = nn.Linear(input_dim, hidden_dim)
        self.c = nn.Parameter(Tensor(1))
        self.d = nn.Parameter(Tensor(1))
        self.c.data.fill_(1.6)
        self.d.data.fill_(torch.log(torch.tensor(19.0)))

        self.act = nn.GELU()
        self.layers = nn.ModuleList([])
        dims = [hidden_dim, hidden_dim]  # 2 layers
        for l in range(len(dims)):
            concat = True #if l==0 else False
            last_dim = hidden_dim if l < len(dims) - 1 else output_dim
            self.layers.append(nn.ModuleList([
                BatchNorm(dims[l]),
                GraphAttention_layer(dims[l], dims[l], heads=heads_num,
                                     concat=concat, dropout=dropout,
                                     attention_type=attention_type,
                                     ),
                GraphAttention_layer(dims[l], dims[l], heads=heads_num,
                                     concat=concat, dropout=dropout,
                                     attention_type=attention_type,
                                     flow='target_to_source',
                                     ),
                nn.Sequential(
                    nn.Linear(dims[l]*2, hidden_dim),
                    nn.GELU(),
                    nn.Dropout(dropout),
                    nn.Linear(hidden_dim, last_dim),
                ),
            ]))
        self.project = nn.Linear(output_dim, output_dim*4)

    def forward(self, data: dict):
        x, edge_index = data['x'], data['edge_index']
        x_auxiliary = torch.sigmoid(self.c * data['node_score_auxiliary'] - self.d)

        x = self.x_input(x)

        att_weights_in, att_weights_out = [], []
        for norm, attn_in, attn_out, ffn in self.layers:
            x = norm(x)
            x_in, att_weights_in_ = attn_in(x, edge_index, x_auxiliary, return_attention_weights=True)
            x_out, att_weights_out_ = attn_out(x, edge_index, x_auxiliary, return_attention_weights=True)
            x = ffn(torch.cat((self.act(x_in), self.act(x_out)), 1))

            att_weights_in.append(att_weights_in_)
            att_weights_out.append(att_weights_out_)
        self.x_embs = x

        # [edge_index, att_in, att_out]
        self.att_weights_first = (att_weights_in[0][0], att_weights_in[0][1], att_weights_out[0][1])
        self.att_weights_second = (att_weights_in[1][0], att_weights_in[1][1], att_weights_out[1][1])

        return self.project(x)


class cl_GRN:
    def __init__(self,
                 hidden_dim: int = 128,
                 output_dim: int = 64,
                 heads_first: int = 4,
                 attention_type: str = 'COS',
                 dropout: float = 0.1,
                 miu: float = 0.5,
                 epochs: int = 500,
                 repeats: int = 5,
                 seed: int = 2022,
                 cuda: int = -1,
                 ):
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.heads_first = heads_first
        self.attention_type = attention_type
        self.dropout = dropout
        self.miu = miu
        self.epochs = epochs
        self.repeats = repeats
        if seed > -1:
            pyg.seed_everything(seed)
            torch.backends.cudnn.deterministic = True
        self.cuda = cuda

        self.idx_GeneName_map = None
        self.att_coefs_ = None
        self.node_embs_ = None
        self.G_predicted = None
        self.att_coef_all = None

    def __get_PYG_data(self, adata: AnnData) -> Data:
        # edge index
        source_nodes = adata.uns['edgelist']['from'].tolist()
        target_nodes = adata.uns['edgelist']['to'].tolist()
        edge_index = torch.tensor([source_nodes, target_nodes], dtype=torch.long)

        # node feature
        x = torch.from_numpy(adata.to_df().T.to_numpy())
        pyg_data = Data(x=x, edge_index=edge_index)

        # node auxiliary score - Differential expression
        if 'node_score_auxiliary' in adata.var:
            pyg_data.node_score_auxiliary = torch.tensor(adata.var['node_score_auxiliary'].to_numpy(),
                                                         dtype=torch.float32).view(-1, 1)
        else:
            print('Warning: Auxiliary gene scores (e.g., differential expression level) are not considered!')
            pyg_data.node_score_auxiliary = torch.ones(x.size(0), dtype=torch.float32).view(-1, 1)

        self.idx_GeneName_map = adata.varm['idx_GeneName_map']

        return pyg_data

    def __corruption(self, data: Data):
        x, edge_index, node_score_auxiliary = data['x'], data['edge_index'], data['node_score_auxiliary']
        data_neg = Data(x=x[torch.randperm(x.size(0))], edge_index=edge_index)
        data_neg.node_score_auxiliary = node_score_auxiliary[torch.randperm(node_score_auxiliary.size(0))]
        return data_neg

    def __summary(self, z, *args, **kwargs):
        #return torch.sigmoid(z.mean(dim=0))
        return torch.sigmoid(torch.cat((3*z.mean(dim=0).unsqueeze(0),
                                        z.max(dim=0)[0].unsqueeze(0),
                                        z.min(dim=0)[0].unsqueeze(0),
                                        2*z.median(dim=0)[0].unsqueeze(0),
                                        ), dim=0))

    def __train(self, data, model, optimizer):
        model.train()
        optimizer.zero_grad()
        pos_z, neg_z, summary = model(data)
        loss = model.loss(pos_z, neg_z, summary)
        loss.backward()
        optimizer.step()
        return float(loss.item())

    def __get_encoder_results(self, data, model):
        model.eval()
        emb_last = model(data)
        return model.x_embs, model.att_weights_first, model.att_weights_second, emb_last

    def run(self, adata: AnnData, showProgressBar: bool=True):
        #device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if self.cuda == -1:
            device = "cpu"
        else:
            device = 'cuda:%s' % self.cuda

        ## Data for pyg input
        data = self.__get_PYG_data(adata).to(device)
        input_dim = data.num_node_features

        ## Run for many times and take the average
        att_weights_all = []
        emb_out_avg = 0
        for iter in range(self.repeats):
            ## Encoder & Model & Optimizer
            encoder = GRN_Encoder(input_dim, self.hidden_dim, self.output_dim, self.heads_first,
                                  dropout=self.dropout, attention_type=self.attention_type).to(device)
            DGI_model = DeepGraphInfomax(hidden_channels=self.output_dim*4,
                                         encoder=encoder,
                                         summary= self.__summary,
                                         corruption=self.__corruption).to(device)
            optimizer = torch.optim.Adam(DGI_model.parameters(), lr=1e-4, weight_decay=5e-4)

            ## Train
            best_encoder = encoder
            min_loss = np.inf
            if showProgressBar:
                with trange(self.epochs, ncols=100) as t:
                    for epoch in t:
                        loss = self.__train(data, DGI_model, optimizer)
                        t.set_description('  Iter: {}/{}'.format(iter + 1, self.repeats))
                        if epoch < self.epochs - 1:
                            t.set_postfix(loss=loss)
                        else:
                            t.set_postfix(loss=loss, min_loss=min_loss)

                        if min_loss > loss:
                            min_loss = loss
                            best_encoder = encoder.state_dict()
            else:
                print('  Iter: {}/{}'.format(iter + 1, self.repeats))
                for epoch in range(self.epochs):
                    loss = self.__train(data, DGI_model, optimizer)
                    if min_loss > loss:
                        min_loss = loss
                        best_encoder = encoder.state_dict()
                print('  Min_Train_loss: {}'.format(min_loss))

            ## Get the result of the best model
            encoder.load_state_dict(best_encoder)
            gene_emb, weights_first, weights_second, emb_last = self.__get_encoder_results(data, encoder)
            gene_emb = gene_emb.cpu().detach().numpy()

            # Use the average of multi-head's attention weights
            weights_first = (weights_first[0].cpu().detach(),
                             torch.cat((weights_first[1].mean(dim=1, keepdim=True),
                                        weights_first[2].mean(dim=1, keepdim=True)), 1).cpu().detach())
            weights_second = (weights_second[0].cpu().detach(),
                              torch.cat((weights_second[1].mean(dim=1, keepdim=True),
                                         weights_second[2].mean(dim=1, keepdim=True)), 1).cpu().detach())

            # Combine the attention coefficients of the first and the second layer
            att_weights = self.miu * weights_first[1] + (1 - self.miu) * weights_second[1]
            att_weights_all.append(att_weights)
            emb_out_avg = emb_out_avg + gene_emb
            if device == 'cuda':
                torch.cuda.empty_cache()

        ## Take the average of multiple runs
        if self.repeats > 1:
            att_weights_all = torch.stack((att_weights_all),0)
        else:
            att_weights_all = att_weights_all[0].unsqueeze(0)
        emb_out_avg = emb_out_avg / self.repeats

        self.edge_index = data.edge_index.cpu()
        self.att_coefs_ = (weights_first[0], att_weights_all) # (edge_index, [n_repeats, edges, 2])
        self.node_embs_ = emb_out_avg

    def get_network(self,
                    keep_self_loops: bool = True,
                    edge_threshold_avgDegree: Optional[int] = 10,
                    edge_threshold_zscore: Optional[float] = None,
                    output_file: Optional[str] = None):
        edge_index_ori = self.edge_index
        edge_index_with_selfloop = self.att_coefs_[0]
        ori_att_coef_all = pd.DataFrame({'from': self.att_coefs_[0][0].numpy().astype(int),
                                         'to': self.att_coefs_[0][1].numpy().astype(int),
                                         'att_coef_in': self.att_coefs_[1].mean(0,keepdim=False)[:,0].numpy(),
                                         'att_coef_out': self.att_coefs_[1].mean(0,keepdim=False)[:,1].numpy()})
        ori_att_coef_all['edge_idx_tmp'] = ori_att_coef_all['from'].astype(str) + "|"\
                                           + ori_att_coef_all['to'].astype(str)

        ## Scaled the attention coefficient for global ranking
        scaled_att_coefficient = []
        g = nx.from_edgelist(edge_index_with_selfloop.numpy().T, create_using=nx.DiGraph)
        for i in range(2):
            # i==0: incoming; i==1: outgoing
            att_weights_i = self.att_coefs_[1][:, :, i] #[num_repeats, num_edges]

            '''
            Edges are selected based on a global weight threshold.
            Attention coefficients are multiplied by the degree of central node.
            '''
            if i == 0:
                d = pd.DataFrame(g.in_degree(), columns=['index', 'degree'])
            else:  # i == 1
                d = pd.DataFrame(g.out_degree(), columns=['index', 'degree'])
            d.index = d['index']
            # att_weights_i * degree_i
            att_weights_i = att_weights_i * np.array(d.loc[edge_index_with_selfloop[1 - i, :].numpy(), 'degree'])
            att_weights_i = att_weights_i.t() #[num_edges, num_repeats]

            if not keep_self_loops:
                edge_index, att_weights_i = remove_self_loops(edge_index_with_selfloop, att_weights_i)
            else:
                # only keep the self loops from the reference network
                data_selfloop_nodes = edge_index_ori[0, edge_index_ori[0] == edge_index_ori[1]]
                selfloop_not_in_data = (edge_index_with_selfloop[0] == edge_index_with_selfloop[1]) & \
                                       ~(edge_index_with_selfloop[0][..., None] == data_selfloop_nodes).any(-1)
                edge_index, att_weights_i = (edge_index_with_selfloop[:, ~selfloop_not_in_data],
                                             att_weights_i[~selfloop_not_in_data])
            # All edge weights without cutoff (num_edge_index, num_repeats)
            scaled_att_coefficient = scaled_att_coefficient + [att_weights_i.clone()]

        ## Output the all edge coefficient of the priori network: 0: incoming; 1: outgoing
        att_weights_combined = (scaled_att_coefficient[0] * 0.5) + (scaled_att_coefficient[1] * 0.5)
        scaled_att_coef_all = pd.DataFrame({'from': edge_index[0].numpy().astype(int),
                                            'to': edge_index[1].numpy().astype(int),
                                            'weights_in': scaled_att_coefficient[0].mean(1,keepdim=False).numpy(),
                                            'weights_out': scaled_att_coefficient[1].mean(1,keepdim=False).numpy(),
                                            'weights': att_weights_combined.mean(1,keepdim=False).numpy(),
                                            'weights_std': att_weights_combined.std(1,keepdim=False).numpy(),
                                            })
        ori_att_coef_all = ori_att_coef_all.loc[ori_att_coef_all['edge_idx_tmp']
                                                    .isin(scaled_att_coef_all['from'].astype(str)
                                                          + "|" + scaled_att_coef_all['to'].astype(str)),].copy()
        self.att_coef_all = pd.merge(ori_att_coef_all[['from', 'to', 'att_coef_in', 'att_coef_out']],
                                     scaled_att_coef_all, on=['from', 'to'], how='inner')
        # Remove the noise edges. Make sure the coefficient of variation (CV) < 0.2.
        CV = scaled_att_coef_all['weights_std'] / (scaled_att_coef_all['weights']+1e-9)
        cv_filter = (CV < 0.2)

        ## Select edges with a cutoff threshold
        att_weights_combined = scaled_att_coef_all['weights']
        if edge_threshold_avgDegree is not None:
            # Select top (N_nodes*edge_threshold_avgDegree) edges
            newEdge_idx = np.argsort(-att_weights_combined)[0:g.number_of_nodes() * edge_threshold_avgDegree]
            newEdge_idx = np.intersect1d(newEdge_idx, np.where(cv_filter)[0])
        else:
            if edge_threshold_zscore is not None:
                # Select top edges with z-score>edge_threshold_zscore
                m, s = att_weights_combined.mean(), att_weights_combined.std()
                edge_threshold_weight = m + (edge_threshold_zscore * s)
                newEdge_idx = np.where((att_weights_combined > edge_threshold_weight) & cv_filter)[0]
            else:  # All edges
                newEdge_idx = list(range(len(att_weights_combined)))
        new_net_df = scaled_att_coef_all.iloc[newEdge_idx,:].copy()

        ## Output the scaled attention coefficient of the predicted network
        ori_att_coef_G = ori_att_coef_all.loc[ori_att_coef_all['edge_idx_tmp']
                                                  .isin(
            new_net_df['from'].astype(str) + "|" + new_net_df['to'].astype(str)),
                                              ['from', 'to', 'att_coef_in', 'att_coef_out']].copy()
        new_net_df = pd.merge(new_net_df, ori_att_coef_G, on=['from', 'to'], how='inner')

        ## To networkx
        G_nx = nx.from_pandas_edgelist(new_net_df, source='from', target='to', edge_attr=True,
                                       create_using=nx.DiGraph)
        num_nodes = len(set(edge_index_with_selfloop[0].numpy()).union(set(edge_index_with_selfloop[1].numpy())))
        largest_components = max(nx.weakly_connected_components(G_nx), key=len)
        if (len(largest_components) / num_nodes) < 0.5:
            print('Warning: the size of maximal connected subgraph is less than half of the input whole graph!')
        G_nx = G_nx.subgraph(largest_components)

        ## Use gene name as index
        mappings = self.idx_GeneName_map.loc[self.idx_GeneName_map['idx'].isin(G_nx.nodes()),
                                             ['idx', 'geneName']]
        mappings = {idx: geneName for (idx, geneName) in np.array(mappings)}
        G_nx = nx.relabel_nodes(G_nx, mappings)
        self.G_predicted = G_nx

        ## Save the predicted network to file
        if isinstance(output_file, str):
            nx.write_edgelist(G_nx, output_file, delimiter=',', data=['weights'])

        return G_nx

    def get_gene_embeddings(self, output_file: Optional[str] = None):
        emb = pd.DataFrame(self.node_embs_, index=self.idx_GeneName_map['geneName'])
        emb = emb.loc[emb.index.isin(self.G_predicted.nodes), :]
        emb.index = emb.index.astype(str)

        if isinstance(output_file, str):
            emb.to_csv(output_file, index_label='geneName')

        return emb

    def get_gene_influence_scores(self):
        gene_influence_scores = pd.DataFrame(np.zeros((len(self.G_predicted.nodes), 2)),
                                             index=sorted(self.G_predicted.nodes),
                                             columns=['out', 'in'])
        for i, v in enumerate(['in', 'out']):
            # The out-degree type of influence is obtained from the incoming network;
            # The in-degree type of influence is obtained from the outgoing network.
            node_att_score = np.sum(nx.to_numpy_matrix(self.G_predicted,
                                                       nodelist=sorted(self.G_predicted.nodes),
                                                       dtype='float32',
                                                       weight='{}_weights'.format(v)),
                                    axis=1 - i)
            gene_influence_scores.iloc[:, i] = np.log1p(node_att_score).flatten().tolist()[0]

        lam = 0.8
        gene_influence_scores['Score'] = lam * gene_influence_scores.loc[:, 'out'] + \
                                         (1-lam) * gene_influence_scores.loc[:, 'in']
        node_weighted_degree = gene_influence_scores.sort_values(by='Score', ascending=False)

        return node_weighted_degree
