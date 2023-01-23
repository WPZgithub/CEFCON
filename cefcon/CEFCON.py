import argparse
from os import fspath
from pathlib import Path

from cellLineage_GRN import cl_GRN
from driver_regulators import driver_regulators, highly_weighted_genes
from utils import *

def main():
    parser = argparse.ArgumentParser(prog='CEFCON', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = add_main_args(parser)
    args = parser.parse_args()

    ## output dir
    p = Path(args.out_dir)
    if not p.exists():
        Path.mkdir(p)

    ## load data
    print('Data loading and preprocessing...')
    data = data_preparation(args.input_expData, args.input_priorNet,
                            genes_DE=args.input_genesDE,
                            TF_list=args.TFs,
                            additional_edges_pct=args.additional_edges_pct)
    if args.input_genesDE is not None:
        genes_DEscore = data.var_names[data.var['node_score_auxiliary']>1]
    else:
        genes_DEscore = None

    ## GRN construction
    print('Constructing cell-lineage-specific GRN...')
    cefcon_GRN_model = cl_GRN(hidden_dim=args.hidden_dim,
                              output_dim=args.output_dim,
                              heads_first=args.heads,
                              attention_type=args.attention,
                              miu=args.miu,
                              epochs=args.epochs,
                              repeats=args.repeats,
                              seed=args.seed,
                              cuda=args.cuda,
                              )
    cefcon_GRN_model.run(data, showProgressBar=True)
    G_predicted = cefcon_GRN_model.get_network(edge_threshold_zscore=None,
                                               edge_threshold_avgDegree=args.edge_threshold_param,
                                               keep_self_loops=~args.remove_self_loops,
                                               output_file=fspath(p/'cl_GRN.csv'))
    node_embeddings = cefcon_GRN_model.get_gene_embeddings(output_file=fspath(p/'gene_embs.csv'))
    gene_influence_scores = cefcon_GRN_model.get_gene_influence_scores()

    ## Driver regulators
    print('Identifying driver regulators...')
    critical_genes, out_critical_genes, in_critical_genes = highly_weighted_genes(gene_influence_scores,
                                                                                  topK=args.topK_drivers)
    cellFate_drivers_set, MDS_driver_set, MFVS_driver_set = driver_regulators(G_predicted,
                                                                              gene_influence_scores,
                                                                              topK=args.topK_drivers,
                                                                              driver_union=True,
                                                                              plot_Venn=False)
    # Driver genes ranking save to file
    drivers_results = gene_influence_scores.loc[gene_influence_scores.index.isin(list(cellFate_drivers_set)), :].copy()
    drivers_results['is_MDS'] = np.isin(drivers_results.index, list(MDS_driver_set))
    drivers_results['is_MFVS'] = np.isin(drivers_results.index, list(MFVS_driver_set))
    if args.TFs is not None:
        TFs = data.var_names[data.var['is_TF'] == 1]
        drivers_results['is_TF'] = np.isin(drivers_results.index, TFs.values)
    drivers_results = drivers_results.sort_values(by='Score', ascending=False)
    drivers_results.to_csv(fspath(p/'driver_regulators.csv'))

    ## RGMs
    print('Identifying regulon-like gene modules...')
    RGMs_results = regulon_activity(data.to_df(), G_predicted,
                                    out_critical_genes.intersection(cellFate_drivers_set),
                                    in_critical_genes.intersection(cellFate_drivers_set),
                                    DEgenes=genes_DEscore,
                                    cell_label=None)
    RGMs = pd.DataFrame([{'Driver_Regulator':r.name, 'Members': list(r.gene2weight.keys())} for r in RGMs_results['regulons']])
    RGMs.to_csv(fspath(p/'RGMs.csv'))

    print('Done. Please check the results in "%s"' % args.out_dir)


def add_main_args(parser: argparse.ArgumentParser):
    # Input data
    input_parser = parser.add_argument_group(title='Input data options')
    input_parser.add_argument('--input_expData', type=str, required=True, metavar='PATH',
                              help='input expression data file')
    input_parser.add_argument('--input_priorNet', type=str, required=True, metavar='PATH',
                              help='input prior network file')
    input_parser.add_argument('--input_genesDE', type=str, default=None, metavar='PATH',
                              help='input differential expression score file')
    input_parser.add_argument('--TFs', type=str, default=None, metavar='PATH',
                              help='input transcriptional factors list')
    input_parser.add_argument('--additional_edges_pct', type=float, default=0.01,
                              help='percentage of additional interactions with highly co-expressions')

    # GRN
    grn_parser = parser.add_argument_group(title='Cell-lineage-specific GRN construction options')
    grn_parser.add_argument('--cuda', type=int, default=0,
                            help="an integer greater than -1 indicates the GPU device number and -1 indicates the CPU device")
    grn_parser.add_argument('--seed', type=int, default=2022,
                            help="random seed (set to -1 means no random seed is assigned)")

    grn_parser.add_argument("--hidden_dim", type=int, default=128,
                            help="hidden dimension of the GNN encoder")
    grn_parser.add_argument("--output_dim", type=int, default=64,
                            help="output dimension of the GNN encoder")
    grn_parser.add_argument("--heads", type=int, default=4,
                            help="number of heads")
    grn_parser.add_argument("--attention", type=str, default='COS', choices=['COS', 'AD', 'SD'],
                            help="type of attention function")
    grn_parser.add_argument('--miu', type=float, default=0.5,
                            help='parameter for considering the importance of attention coefficients of the first GNN layer')
    grn_parser.add_argument('--epochs', type=int, default=350,
                            help='number of epochs')
    grn_parser.add_argument('--repeats', type=int, default=5,
                            help='number of repeats')

    grn_parser.add_argument("--edge_threshold_param", type=int, default=8,
                            help="threshold for selecting top-weighted edges")
    grn_parser.add_argument("--remove_self_loops", action="store_true",
                            help="remove self loops")

    # Driver regulators
    driver_parser = parser.add_argument_group(title='Driver regulator identification options')
    driver_parser.add_argument('--topK_drivers', type=int, default=50,
                               help="number of candidate drivers genes according to the influence score")

    # Output dir
    parser.add_argument("--out_dir", type=str, required=True, default='./output',
                        help="results output path")

    return parser

if __name__ == "__main__":
    main()
