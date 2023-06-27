"""

Algorithmic description
-----------------------

The procedure is explained in detail on the following paper, whose datasets can
be found in the corresponding repository:

López-Rozo, N.; Ramirez-Castrillon, M.; Romero, M.; Finke, J.; Rocha, C.
*Gene Expression Datasets for Two Versions of the Saccharum spontaneum AP85-441 Genome*.
Data 2023, 8, 1. https://doi.org/10.3390/data8010001

Based on the output of BLASTN, the associations among the alleles in v2018 and
v2019 are found to have repetitions. In the case of the mappings between v2018
to v2019, a CDS in the source could be associated with several CDS in the
target. To generate a reasonable coverage, both mappings are combined by
modeling the problem as a graph flow optimization problem [20] with multiple
sources (v2018 alleles) and multiple targets (v2019 alleles).

A min-cost max-flow problem requires to compute a graph-matching (i.e., match
at the level of nodes/vertices) with maximal cardinality (i.e., maximal number
of connections), thus ensuring a maximal covering of the source-target
associations. If more than one maximal matching is possible, then the cost of
producing that maximal flow is to be minimized. In this case, identity scores
can be considered to identify the matching with the greatest sum of identity
scores, while still ensuring that a v2018 allele expression is used at most
once. Since the algorithm implemented in networkx minimizes cost, the
artificial cost fed to the min-cost max-flow algorithm is *pident* (i.e.,
percent identity) on each possible association between the two versions of the
alleles.

|

**API Reference**
-----------------
"""

import pickle as pk
import time

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_expression(file, sep=None):
    """
    This function reads the expression file and returns a list of names and a
    matrix of expression values.

    :param file: The path to the expression file.
    :type file: str
    :param sep: The separator used in the file.
    :type sep: str

    :return names: The list of names.
    :rtype names: list
    :return mat: The matrix of expression values.
    :rtype mat: list[numpy.array]
    :return headers: The headers of the table ("gene"/"allele" + accession names).
    :rtype headers: str
    """
    names = []
    mat = []
    with open(file, "r") as f:
        headers = f.readline()
        for line in f.readlines():
            name, tmp = line.strip().split(sep, 1)
            names.append(name)
            mat.append(np.array([float(x) for x in tmp.strip().split(sep)]))
    return names, mat, headers


def create_bipartite_graph(df_map, query_col='qseqid', subject_col='sseqid',
                           weight_col='pident', base_graph=None):
    """
    This function creates a bipartite graph from the mapping dataframe, which was
    extracted from the BLASTN output. **Note**: The weight of the edges is
    multiplied by -1000 and truncated to be used as a cost in the min-cost
    max-flow algorithm.

    :param df_map: The mapping dataframe.
    :type df_map: pandas.DataFrame
    :param query_col: The name of the column containing the query sequences.
    :type query_col: str
    :param subject_col: The name of the column containing the subject sequences.
    :type subject_col: str
    :param weight_col: The name of the column containing the weight of the edges.
    :type weight_col: str
    :param base_graph: The base graph to be used. if set to None (default), a new \
    graph is created.
    :type base_graph: networkx.Graph

    :return bipartite: The bipartite graph.
    :rtype bipartite: networkx.Graph
    """
    if base_graph is None:
        bipartite = nx.Graph()
    else:
        bipartite = base_graph
    bipartite.add_nodes_from(df_map[query_col].tolist(), subset=query_col)
    bipartite.add_nodes_from(df_map[subject_col].tolist(), subset=subject_col)
    for index, row in df_map.iterrows():
        if not bipartite.has_edge(row[query_col], row[subject_col]):
            bipartite.add_edge(row[query_col], row[subject_col], capacity=1,
                               weight=int(row[weight_col] * -1000))
        elif bipartite[row[query_col]][row[subject_col]]['weight'] > int(row[weight_col] * -1000):
            bipartite[row[query_col]][row[subject_col]]['weight'] = int(row[weight_col] * -1000)
            # print(bipartite[row[query_col]][row[subject_col]])
    return bipartite


def draw_multigraph(graph, subset_feat, x_offset=1, y_offset=1,
                    subset_color=("green", "blue", "purple", "red", "orange", "yellow"),
                    savefig=False, filename="graph.pdf", node_size=100):
    """
    This function draws a multipartite graph. **Note**: The graph is drawn in
    the order of the subsets in the list of subset features. The nodes are

    :param graph: The graph to be drawn.
    :type graph: networkx.Graph
    :param subset_feat: The subset feature names. The order of the names is the order \
    the subsets will be drawn.
    :param x_offset: The x offset for the subsets and the drawn nodes.
    :type x_offset: int
    :param y_offset: The y offset for the subsets and the drawn nodes.
    :type y_offset: int
    :param subset_color: The colors of the subsets.
    :type subset_color: list[str]
    :param savefig: Whether to save the figure or not.
    :type savefig: bool
    :param filename: The name of the file to save the figure.
    :type filename: str
    :param node_size: The size of the nodes.
    :type node_size: int
    """
    # Extract the nodes of the graph, according to the subsets
    subset_nodes = []
    subset_sizes = []
    colors = [subset_color[i % 6] for i in range(len(subset_feat))]
    for subset in subset_feat:
        subset_nodes.append([n for n, d in graph.nodes(data=True) if d['subset'] == subset])
        subset_sizes.append(len([n for n, d in graph.nodes(data=True) if d['subset'] == subset]))

    # compute the node positions, according to the subsets. Each layer will be drawn
    # in a different x coordinate, and the nodes will be centered according to the x-axis.
    node_positions = []
    for i, size in enumerate(subset_sizes):
        node_positions.extend([(i * x_offset, y_offset * (-j * 2 + size - 1)/2) for j in range(size)])

    # compute the accumulated sum of the subset sizes
    s = [0] + list(np.cumsum(subset_sizes))

    # gather a map from the nodes to their positions in the node_positions list
    node_pos_map = dict()
    for i, subset in enumerate(subset_nodes):
        for j, node in enumerate(subset):
            node_pos_map[node] = s[i] + j

    # draw the nodes using matplotlib.pyplot
    ct = 0
    for i, subset in enumerate(subset_nodes):
        for j in subset:
            color = colors[i]
            pos = node_positions[ct]
            plt.scatter(*pos, s=node_size, c=color, edgecolor='black', zorder=1)
            ct += 1

    # draw the edges using matplotlib.pyplot, locating them behind the nodes
    for u, v in graph.edges():
        pos_u = node_positions[node_pos_map[u]]
        pos_v = node_positions[node_pos_map[v]]
        pos = zip(pos_u, pos_v)
        plt.plot(*pos, c='black', alpha=0.5, zorder=0)

    plt.axis('off')
    plt.tight_layout()
    if savefig:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()


def allele_transform(mapping_1, mapping_2, expression_1, result_folder="./",
                     expression_2="expression2.csv", prefix_restricted="restricted",
                     blastn_titles="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,"
                                   "qend,sstart,send,evalue,bitscore",
                     right_subset="qseqid", verbose=False):
    """
    This function transforms the expression values of the alleles in v2018 to
    the alleles in v2019. The transformation is done by solving a min-cost
    max-flow problem with multiple sources (v2018 alleles) and multiple targets
    (v2019 alleles).

    Keep in mind that the function stores additional files in the result folder:

    - *expression_2*: The mapped expression file for v2019, as CSV format.

    - *prefix_restricted.pk*: The gene assignment, as a dictionary.

    - *prefix_restricted.txt*: The gene assignment, as a table.

    :param mapping_1: The path to the mapping file from v2018 to v2019.
    :type mapping_1: str
    :param mapping_2: The path to the mapping file from v2019 to v2018.
    :type mapping_2: str
    :param expression_1: The path to the expression file of v2018.
    :type expression_1: str
    :param result_folder: The path to the folder where the results will be saved.
    :type result_folder: str
    :param expression_2: The name of the expression file of v2019.
    :type expression_2: str
    :param prefix_restricted: The prefix of the restricted alleles.
    :type prefix_restricted: str
    :param blastn_titles: The titles of the BLASTN output.
    :type blastn_titles: str
    :param right_subset: The subset identifier of the left nodes of the bipartite graph.
    :type right_subset: str
    :param verbose: If set to True, the function prints the time at the beginning \
    and the end of the execution.
    :type verbose: bool

    :return bipartite: The bipartite graph used for computing the min-cost max-flow.
    :rtype bipartite: networkx.Graph
    """

    if verbose:
        print(time.strftime("%a, %d %b %Y %H:%M:%S"))
    names, exp, headers = read_expression(expression_1, ",")
    titles = blastn_titles.strip().split(",")
    df_map = pd.read_csv(mapping_1, sep=None, index_col=None,
                         names=titles, engine='python')
    # take duplicated queries
    # df_dup = df_map[df_map.duplicated('sseqid')]
    # rep_CDS = df_dup['sseqid'].unique()
    # tmp = df_map[df_map['sseqid'].isin(rep_CDS)].sort_values(by=['sseqid'])
    # if verbose:
    #     print(tmp)

    # Bipartite graph creation: Phase 1
    bipartite = create_bipartite_graph(df_map, query_col='qseqid', subject_col='sseqid',
                                       weight_col='pident')
    if verbose:
        print(f"Phase 1: (V,E) = ({bipartite.number_of_nodes()}, {bipartite.number_of_edges()})")

    # now load the opposite direction, but change the first two column names
    # ["sseqid", "qseqid", "pident", ...]
    df_map = pd.read_csv(mapping_2, sep=None, index_col=None,
                         names=[titles[1]] + [titles[0]] + titles[2:], engine='python')

    # Note: edge direction is opposite to the previous function call
    bipartite = create_bipartite_graph(df_map, query_col='qseqid', subject_col='sseqid',
                                       weight_col='pident', base_graph=bipartite)
    if verbose:
        print(f"Phase 2: (V,E) = ({bipartite.number_of_nodes()}, {bipartite.number_of_edges()})")
    # "qseqid", "sseqid", "pident", ...
    df_map = pd.read_csv(mapping_1, sep=None, index_col=None, names=titles, engine='python')
    # MAXIMUM MATCHING OF EACH COMPONENT OF BIPARTITE
    # (computed only when the number of nodes to the right is >1)
    restricted = dict()  # dictionary of restricted edges (used in the max flow)
    l_nodes = set([n for n, d in bipartite.nodes(data=True) if d["subset"] == titles[0]])
    r_nodes = set(bipartite.nodes()) - l_nodes

    # Configuring source and target for graph flow problem
    bipartite.add_node('SS', subset='SS')
    bipartite.add_node('ST', subset='ST')
    for node in l_nodes:
        bipartite.add_edge('SS', node, capacity=1, weigth=0)  # zero weight
    for node in r_nodes:
        bipartite.add_edge(node, 'ST', weight=0)  # infinite capacity, zero weight
    maxflow = nx.max_flow_min_cost(bipartite, 'SS', 'ST')

    if verbose:
        print(time.strftime("%a, %d %b %Y %H:%M:%S"), "Sum of PIDENT values:", -0.001*nx.cost_of_flow(bipartite, maxflow))
    ct = [0, 0, 0]
    for i, v in enumerate(bipartite['SS']):
        ct[maxflow['SS'][v]] += 1
    if verbose:
        print("Flow [0, 1, -1] =", ct)
    ct = 0  # counter of the number of assigned genes
    for i, v in enumerate(bipartite['SS']):
        # if maxflow['SS'][v]['flow'] == 1:
        if maxflow['SS'][v] == 0:
            continue
        ct += 1
        for u in maxflow[v]:
            if u == 'SS':
                continue
            # if maxflow[v][u]['flow'] == 1:
            if maxflow[v][u] == 1 and maxflow[u][v] == 0:
                # Flow from v to u goes from left to right only
                if v in restricted:
                    if verbose:
                        print(f"{i}: {v} appears assigned twice! {u} and {restricted[v]}")
                        for x in bipartite[v]:
                            print(f"{v} -> {x}: {maxflow[v][x]}\t {x} -> {v}: {maxflow[x][v]}")
                            for y in bipartite[x]:
                                if y in l_nodes:
                                    continue
                                print(f"    {x} -> {y}: {maxflow[x][y]}\t {y} -> {x}: {maxflow[y][x]}")
                    raise Exception(f"{v} appears assigned twice! {u} and {restricted[v]}")
                else:
                    restricted[v] = u
    # extracting the names of the genes as all nodes adjacent to the 'ST' node
    new_names = [x for x in bipartite['ST']]
    dict_old = {v: u for u, v in enumerate(names)}
    dict_new = {v: u for u, v in enumerate(new_names)}
    mat = np.zeros((len(new_names), len(exp[0])))
    # print(dict_old)
    # print(dict_new)
    if verbose:
        print("Number of resulting transcripts:", ct)
        print("Number of genes:", len(new_names))
    for key in restricted:
        value = restricted[key]
        # print(key, value)
        mat[dict_new[value]] += exp[dict_old[key]]

    # Saving the resulting DataFrame
    with open(result_folder + expression_2, "w") as f:
        f.write(headers)
        for i, name in enumerate(new_names):
            f.write(name)
            for x in mat[i]:
                f.write("," + str(float(x)))
            f.write("\n")
    pk.dump(restricted, open(f"{result_folder}{prefix_restricted}.pk", "wb"))
    with open(f"{result_folder}{prefix_restricted}.txt", "w") as f:
        f.write("key,value\n")
        for key in restricted:
            f.write(f"{key},{restricted[key]}\n")
    if verbose:
        print(time.strftime("%a, %d %b %Y %H:%M:%S"), "DONE!")

    return bipartite


def _test1():
    """
    Testing the allele_transform function using extracted data from the data folder.
    It plots the resulting graph.
    """
    file_mapping_1 = "blast_1.txt"
    file_mapping_2 = "blast_2.txt"
    file_expression_in = "expression.csv"
    file_expression_out = "expression_2.csv"
    folder_in = "data/"
    folder_out = "test_results/"
    prefix_restricted = "restricted"
    g = allele_transform(folder_in + file_mapping_1, folder_in + file_mapping_2,
                         folder_in + file_expression_in, folder_out,
                         file_expression_out, prefix_restricted, verbose=True)
    print("Results from the allele transform algorithm (first 10 lines):")
    with open(f"{folder_out}{prefix_restricted}.txt", "r") as f:
        for line in f.readlines()[:10]:
            print(line.strip())

    # Plotting the resulting graph
    draw_multigraph(g, ["SS", "qseqid", "sseqid", "ST"], x_offset=10,
                    savefig=True, filename="test_results/graph.pdf")

    # pos = nx.spring_layout(g, seed=42)  # Define the layout of the graph
    pos = nx.bipartite_layout(g, list(g['SS']) + ['SS'])  # Define the layout of the graph
    pos['SS'][0] -= 1
    pos['ST'][0] += 1
    edge_labels = nx.get_edge_attributes(g, "capacity")  # Get the edge capacities as labels

    # Draw nodes and edges
    nx.draw_networkx_nodes(g, pos, node_color="lightblue", node_size=10)
    nx.draw_networkx_edges(g, pos, edge_color="gray")
    # nx.draw_networkx_labels(g, pos)
    nx.draw_networkx_edge_labels(g, pos, edge_labels=edge_labels)

    # Display the graph
    plt.axis("off")
    plt.show()


if __name__ == '__main__':
    _test1()
