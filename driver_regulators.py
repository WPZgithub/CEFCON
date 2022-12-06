import networkx as nx
import numpy as np
import pandas as pd
import cvxpy as cvx

def _root_variables(directed_graph):
    return set([n for n in directed_graph.nodes() if (directed_graph.in_degree(n) == 0)
                or (list(directed_graph.predecessors(n)) == [n])])


def _end_variables(directed_graph):
    return set([n for n in directed_graph.nodes() if (directed_graph.out_degree(n) == 0)
                or (list(directed_graph.successors(n)) == [n])])


"""
MDS control for directed network (Graph reduction + ILP)
Ref: [1] Dominating scale-free networks with variable scaling exponent: heterogeneous networks are not difficult to control
     New Journal of Physics, 2012
     [2] Critical controllability analysis of directed biological networks using efficient graph reduction.
     2017, Scientific Reports
A set S⊆V of nodes in a graph G=(V,E) is a dominating set if every node v∈V is either an element of S 
or adjacent to an element of S.
Use efficient graph reduction and then solve the ILP problem using gurobi. 
"""
def _MDS_graph_reduction(directed_graph):
    critical_nodes = set()
    redundant_nodes = set()

    # Critical nodes cond. 1: Source nodes are critical (driver) nodes
    critical_nodes.update(_root_variables(directed_graph))
    remain_nodes = set(directed_graph.nodes()) - critical_nodes

    noChange = False
    while not noChange:
        noChange = True

        # Critical nodes cond. 2: A node with at least two directed edges to nodes with outdegree 0 and indegree 1
        in1out0_nodes = set([n for n in directed_graph.nodes() if (directed_graph.in_degree(n) == 1
                                                                  and directed_graph.out_degree(n) == 0)])
        add_critical = set([n for n in list(remain_nodes-in1out0_nodes)
                            if (len(set(directed_graph.successors(n)).intersection(in1out0_nodes))>1)])
        if len(add_critical)==0:
            noChange *= True
        else:
            critical_nodes.update(add_critical)
            remain_nodes = remain_nodes - add_critical
            remove_edges = [(i,n) for n in list(add_critical) for i in directed_graph.predecessors(n)]
            directed_graph.remove_edges_from(remove_edges)
            noChange *= False

        # Redundant nodes cond.: a node with outdegree 0 and has an incoming link from a critical node
        add_redundant = set([n for n in list(remain_nodes)
                             if (len(set(directed_graph.predecessors(n)).intersection(critical_nodes))>0
                                 and (directed_graph.out_degree(n)==0))])
        if len(add_redundant)==0:
            noChange *= True
        else:
            redundant_nodes.update(add_redundant)
            remain_nodes = remain_nodes - add_redundant
            directed_graph.remove_nodes_from(add_redundant)
            noChange *= False

    return directed_graph, critical_nodes, redundant_nodes

def MDScontrol(directed_graph, solver='GUROBI'):
    print('  Solving MDS problem...')
    directed_graph.remove_edges_from(nx.selfloop_edges(directed_graph))
    reduced_graph = nx.DiGraph(directed_graph)
    intermittent_nodes = set()
    MDS_driver_set = set()

    ## Graph reduction
    reduced_graph, critical_nodes, redundant_nodes = _MDS_graph_reduction(reduced_graph)
    print('    {} critical nodes are found.'.format(len(critical_nodes)))

    ## ILP is used to find an MDS in the reduced graph
    reduced_graph.remove_nodes_from(list(nx.isolates(reduced_graph)))
    print('    {} nodes left after graph reduction operation.'.format(reduced_graph.number_of_nodes()))
    if reduced_graph.number_of_nodes() == 0:
        print('{} MDS driver nodes are found.'.format(len(critical_nodes)))
    else:
        print('    Solving the Integer Linear Programming problem on the reduced graph...')
        A = nx.to_numpy_matrix(reduced_graph)
        A = A + np.diag(np.ones(A.shape[0]))
        # Define the optimization variables
        x = cvx.Variable(reduced_graph.number_of_nodes(), boolean=True)
        # Define the constraints
        constraints = [A@x >= np.ones(reduced_graph.number_of_nodes())]
        # Define the optimization problem
        obj = cvx.Minimize(cvx.sum(x))
        # Solve the problem
        prob = cvx.Problem(obj, constraints)

        if solver=='GUROBI':
            # Solve with GUROBI.
            import gurobipy as gp
            m = gp.Model()
            print('      Solving by GUROBI...(', end='')
            prob.solve(solver=cvx.GUROBI, verbose=False)
            print('optimal value with GUROBI:{},'.format(prob.value), end='  ')
        #elif solver == 'XPRESS':
        #        prob.solve(solver=cvx.XPRESS, verbose=False)
        #        print("optimal value with XPRESS:", prob.value)
        else:
            # Solve with ECOS_BB
            print('      Inaccurate solver is selected! Now, solving by SCIP...(', end='')
            prob.solve(solver=cvx.SCIP, verbose=False)
            print('optimal value with SCIP:{},'.format(prob.value), end='  ')
        print('status:{})'.format(prob.status))

        ## Set remain nodes that belongs to the MDS as critical nodes
        nodes_idx_map = dict(zip(range(reduced_graph.number_of_nodes()), reduced_graph.nodes()))
        mds_nodes = set([v for k,v in nodes_idx_map.items() if x.value[k]==1])
        MDS_driver_set = critical_nodes.union(mds_nodes)
        intermittent_nodes = set(reduced_graph.nodes()) - MDS_driver_set
        print('  {} MDS driver genes are found.'.format(len(MDS_driver_set)))

    return MDS_driver_set, intermittent_nodes

"""
FVS control for directed network (Graph reduction + ILP)
Ref: [1] Structure-based control of complex networks with nonlinear dynamics.
     Proceedings of the National Academy of Sciences, 2017
     [2] On computing the minimum feedback vertex set of a directed graph by contraction operations.
     IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 2000
     [ILP] An exact algorithm for selecting partial scan flip-flops. 1994
Use efficient graph reduction and then solve the ILP problem using gurobi. 
"""
def _in0out0(directed_graph):
    remove_set = set([n for n in directed_graph.nodes() if (directed_graph.in_degree(n) == 0
                                                            or directed_graph.out_degree(n) == 0)])
    directed_graph.remove_nodes_from(remove_set)
    return directed_graph, len(remove_set) == 0

def _selfloop(directed_graph, S):
    remove_set = list(nx.nodes_with_selfloops(directed_graph))
    S = S.union(set(remove_set))
    directed_graph.remove_nodes_from(remove_set)
    return directed_graph, S, len(remove_set) == 0

def _in1(directed_graph):
    remove_set = set([n for n in directed_graph.nodes() if (directed_graph.in_degree(n) == 1)])
    temp_set = set()
    for u in remove_set:
        v = list(directed_graph.predecessors(u))[0]
        if v != u:
            directed_graph.add_edges_from([(v, w) for w in directed_graph.successors(u)])
            directed_graph.remove_node(u)
        else:
            temp_set.update(u)
    remove_set = remove_set - temp_set
    return directed_graph, len(remove_set) == 0

def _out1(directed_graph):
    remove_set = set([n for n in directed_graph.nodes() if (directed_graph.out_degree(n) == 1)])
    temp_set = set()
    for u in remove_set:
        v = list(directed_graph.successors(u))[0]
        if u!= v:
            directed_graph.add_edges_from([(w, v) for w in directed_graph.predecessors(u)])
            directed_graph.remove_node(u)
        else:
            temp_set.update(u)
    remove_set = remove_set - temp_set
    return directed_graph, len(remove_set) == 0

def _PIE(directed_graph):
    G_SCCs = [c for c in nx.strongly_connected_components(directed_graph)]
    edges_SCCs = set()
    for g in G_SCCs:
        edges_SCCs.update(set(directed_graph.subgraph(g).edges()))
    remove_set = set(directed_graph.edges()) - edges_SCCs
    directed_graph.remove_edges_from(list(remove_set))
    return directed_graph, len(remove_set) == 0

def _CORE(directed_graph, nodes_importance, S):
    remove_set = set()
    pie = [e for e in directed_graph.edges() if (directed_graph.has_edge(e[1],e[0]) and (e[1]!=e[0]))]
    if len(pie) > 0:
        x,y = zip(*pie)
        nodes_pie = set(x).union(set(y))
        piv = [n for n in list(nodes_pie)
               if (set(directed_graph.predecessors(n))==set(directed_graph.successors(n)))]
        #sort PIV in the ascending order according to the degree and out_degree_importance
        piv_df = pd.DataFrame(directed_graph.degree(piv), index=piv, columns=['id','degree'])
        piv_df['importance'] = piv_df.index.map(nodes_importance)
        piv_df.sort_values(by=['importance','degree'], ascending=[True,True], inplace=True)
        #piv_df.sort_values(by='degree', ascending=True, inplace=True)
        piv_df['valid'] = True
        for v in piv_df.index:
            if piv_df.loc[v,'valid']:
                # d-clique
                d_graph = directed_graph.subgraph([v] + list(directed_graph.neighbors(v)))
                is_dclique = True
                for d_i in d_graph:
                    if len(list(nx.bfs_edges(d_graph, d_i, depth_limit=1)))<(len(d_graph)-1):
                        is_dclique = False

                if is_dclique:
                    S = S.union(set(d_graph.nodes()) - set(v)) # v is CORE
                    remove_set = remove_set.union(set(d_graph.nodes()))
                    piv_df.loc[piv_df.index.isin(d_graph.nodes()), 'valid'] = False
                else:
                    piv_df.loc[v, 'valid'] = False
    directed_graph.remove_nodes_from(remove_set)
    return directed_graph, S, len(remove_set) == 0

def _DOME(directed_graph):
    pie = [e for e in directed_graph.edges() if (directed_graph.has_edge(e[1],e[0]) and (e[1]!=e[0]))]
    auxilliary_graph = nx.DiGraph(directed_graph)
    auxilliary_graph.remove_edges_from(pie)
    remove_set = [e for e in auxilliary_graph.edges() if e[0]!=e[1] and
                  (set(auxilliary_graph.predecessors(e[0])).issubset(set(auxilliary_graph.predecessors(e[1])))
                  or
                  set(auxilliary_graph.successors(e[1])).issubset(set(auxilliary_graph.successors(e[0]))))]
    directed_graph.remove_edges_from(remove_set)
    return directed_graph, len(remove_set) == 0

def _MFVS_graph_reduction(directed_graph, nodes_importance, S):
    all_done = False
    while not all_done:
        all_done = True

        directed_graph, no_change = _in0out0(directed_graph)
        all_done *= no_change

        directed_graph, S, no_change = _selfloop(directed_graph, S)
        all_done *= no_change

        directed_graph, no_change = _out1(directed_graph)
        all_done *= no_change

        directed_graph, no_change = _in1(directed_graph)
        all_done *= no_change

        directed_graph, no_change = _PIE(directed_graph)
        all_done *= no_change

        directed_graph, S, no_change = _CORE(directed_graph, nodes_importance, S)
        all_done *= no_change

        directed_graph, no_change = _DOME(directed_graph)
        all_done *= no_change

    return directed_graph,S

def MFVScontrol(directed_graph, nodes_importance, solver='GUROBI'):
    print('  Solving MFVS problem...')
    # Source nodes are critical (steering) nodes
    critical_nodes = set()
    source_nodes = _root_variables(directed_graph)
    critical_nodes.update(source_nodes)
    reduced_graph = nx.DiGraph(directed_graph)

    while True:
        ## Graph reduction
        reduced_graph, critical_nodes = _MFVS_graph_reduction(reduced_graph, nodes_importance, critical_nodes)
        reduced_graph_SCCs = [c for c in nx.strongly_connected_components(reduced_graph)]
        if len(reduced_graph_SCCs) > 1:
            reduced_graph_new = nx.DiGraph()
            for g in reduced_graph_SCCs:
                if len(g) > 1:
                    g = nx.DiGraph(reduced_graph.subgraph(g))
                    g, critical_nodes = _MFVS_graph_reduction(g, nodes_importance, critical_nodes)
                    reduced_graph_new.update(g)
            reduced_graph = nx.DiGraph(reduced_graph_new)

        ## Some tricks (Let the node with the largest degree be a critical node, and do graph reduction again)
        ## In order to facilitate the ILP solution
        if reduced_graph.number_of_nodes() > 150:
            d_seq = pd.merge(pd.DataFrame(reduced_graph.in_degree(),columns=['Node','in_Degree']),
                             pd.DataFrame(reduced_graph.out_degree(),columns=['Node','out_Degree']),
                             on='Node')
            d_seq['dot_Degree'] = d_seq['in_Degree'] * d_seq['out_Degree']
            node_max = d_seq['Node'][d_seq['dot_Degree'].idxmax()]
            critical_nodes.add(node_max)
            reduced_graph.remove_node(node_max)
        else:
            break
    print('    {} critical nodes are found.'.format(len(critical_nodes)))

    ## ILP is used to find an MFVS in the reduced graph
    reduced_graph.remove_nodes_from(list(nx.isolates(reduced_graph)))
    print('    {} nodes left after graph reduction operation.'.format(reduced_graph.number_of_nodes()))
    if reduced_graph.number_of_nodes() == 0:
        print('{} MFVS driver genes are found.'.format(len(critical_nodes)))
    else:
        print('    Solving the Integer Linear Programming problem on the reduced graph...')
        nodes_idx_map = dict(zip(reduced_graph.nodes(), range(len(reduced_graph.nodes()))))
        # Define the optimization variables
        n = reduced_graph.number_of_nodes()
        x = cvx.Variable(n, boolean=True)
        # Define the constraints
        constraints = []
        w = cvx.Variable(n, integer=True)
        for e in reduced_graph.edges():
            constraints += [w[nodes_idx_map[e[0]]] - w[nodes_idx_map[e[1]]] + n*x[nodes_idx_map[e[0]]] >=1]
        constraints += [0 <= w, w <= (n-1)]
        # Define the optimization problem
        obj = cvx.Minimize(cvx.sum(x))
        # Solve the problem
        prob = cvx.Problem(obj, constraints)

        if solver=='GUROBI':
            # Solve with GUROBI.
            import gurobipy as gp
            m = gp.Model()
            print('      Solving by GUROBI...(', end='')
            prob.solve(solver=cvx.GUROBI, verbose=False)
            print('optimal value with GUROBI:{},'.format(prob.value), end='  ')
        #elif solver=='XPRESS':
        #    prob.solve(solver=cvx.XPRESS, verbose=False)
        #    print("optimal value with XPRESS:", prob.value)
        else:
            # Solve with ECOS_BB
            print('      Inaccurate solver is selected. Now, solving by SCIP...(', end='')
            prob.solve(solver=cvx.SCIP, verbose=False)
            print('optimal value with SCIP:{},'.format(prob.value), end='  ')
        print('status:{})'.format(prob.status))

        ## Set remain nodes that belongs to the MDS as critical nodes
        mFVS_remain = set([k for k,v in nodes_idx_map.items() if x.value[v]==1])
        critical_nodes = critical_nodes.union(mFVS_remain)
        print('  {} MFVS driver nodes are found.'.format(len(critical_nodes)))

    return critical_nodes, source_nodes


def highly_weighted_genes(gene_influence_scores, topK_drivers=50):
    v_out =  gene_influence_scores.sort_values(by='out', ascending=False)
    v_out = v_out.iloc[0:topK_drivers,[0]]
    out_critical_genes = set(v_out[v_out>0].index)

    v_in = gene_influence_scores.sort_values(by='in', ascending=False)
    v_in = v_in.iloc[0:topK_drivers,[1]]
    in_critical_genes = set(v_in[v_in>0].index)

    critical_genes = out_critical_genes.union(in_critical_genes)
    return critical_genes, out_critical_genes, in_critical_genes


def driver_regulators(G, gene_influence_scores,
                      topK_threshold=50, driver_union=True,
                      solver='GUROBI', plot_Venn=False):
    DFVS_driver_set, source_nodes = MFVScontrol(G, gene_influence_scores.loc[:, 'out'], solver=solver)
    DFVS_driver_set = DFVS_driver_set.union(source_nodes)
    MDS_driver_set, MDS_intermittent_nodes = MDScontrol(G, solver=solver)
    if driver_union:
        driver_set = MDS_driver_set.union(DFVS_driver_set)
    else:
        driver_set = MDS_driver_set.intersection(DFVS_driver_set)

    critical_genes, _, _ = highly_weighted_genes(gene_influence_scores, topK_threshold)
    CEFCON_drivers = driver_set.intersection(critical_genes)

    # Plot Venn figure
    if plot_Venn:
        from matplotlib import pyplot as plt
        from matplotlib_venn import venn3
        plt.subplots(figsize=(12, 9), dpi=300)
        out = venn3(subsets=[MDS_driver_set, DFVS_driver_set, critical_genes],
                    set_labels=('MDS driver genes({})'.format(len(MDS_driver_set)),
                                'MFVS driver genes({})'.format(len(DFVS_driver_set)),
                                'Top ranked genes({})'.format(len(critical_genes))),
                     set_colors=('#0076FF','#D74715','#009000'),
                     alpha=0.45)
        for text in out.set_labels:
            if text is not None:
                text.set_fontsize(30)
        for text in out.subset_labels:
            if text is not None:
                text.set_fontsize(32)
        plt.tight_layout()
        #plt.savefig(Venn_fig_file+'.svg', dpi=600, format='svg', bbox_inches='tight')
        plt.show()

    return CEFCON_drivers, MDS_driver_set, DFVS_driver_set
