import pickle
import random
import sys

import networkx as nx
from matplotlib import pyplot as plt
from ppidm_validation.rewire_gold_standard import sorted_node_degrees

random.seed(42)


def rewire_graph(ext_graph: nx.Graph, file_path=None):
    sorted_ext_graph = sorted_node_degrees(list(ext_graph.edges))
    nodes = list(sorted_ext_graph.keys())
    degrees = list(sorted_ext_graph.values())

    random_ext_graph = nx.expected_degree_graph(degrees, seed=42)

    node_mapping = {i: j for i, j in zip(range(len(nodes)), nodes)}
    random_ext_graph: nx.Graph = nx.relabel_nodes(random_ext_graph, node_mapping)
    print("Rewired graph:", random_ext_graph)
    if file_path is not None:
        pickle.dump(random_ext_graph, open(file_path, 'wb'))

    random_degrees = {}
    random_interactions = list(random_ext_graph.edges)
    for key in sorted_ext_graph.keys():
        random_degrees[key] = len(list(random_ext_graph.neighbors(key)))

    print("Network Density of random:", nx.density(random_ext_graph))
    print("Network Density of real  :", nx.density(ext_graph))
    # practially the same

    # print("Average Clustering Coefficient of random:", nx.average_clustering(random_ext_graph))
    # print("Average Clustering Coefficient of real  :", nx.average_clustering(ext_graph))
    # random: 0.05019420725696134
    # real: 0.38526512773920396
    return random_ext_graph, degrees, random_degrees


def visualize(degrees, random_degrees, file_path=None, title='Expected node degree vs. true node degree'):
    # Visualising
    plt.style.use('ggplot')
    x = degrees
    y = list(random_degrees.values())

    plt.scatter(x, y, color='#1F77B4')
    # Adding labels to the plot
    plt.xlabel('expected degree')
    plt.ylabel('true degree')
    plt.plot([0, x[0]], [0, x[0]], 'r--')
    plt.title(title)

    # Displaying the plot
    if file_path is not None:
        plt.savefig(file_path)
    plt.show()


def reassign_ddis(ppi_ddi_graph: nx.Graph, random_graph: nx.Graph):

    # find all relevant domains for any given gene/protein
    ddi_pairs = dict()
    for node in ppi_ddi_graph.nodes:
        prot, domain = node.split("/")
        if prot not in ddi_pairs:
            ddi_pairs[prot] = set()
        ddi_pairs[prot].add(domain)

    # assign each ppi two random domains
    domainG_random = set()
    skips = 0
    for prot_a, prot_b in random_graph.edges:
        try:
            prot_a_domain = random.choice(list(ddi_pairs.get(prot_a)))
            prot_b_domain = random.choice(list(ddi_pairs.get(prot_b)))
        except TypeError:
            skips += 1
            continue
        domainG_random.add((f"{prot_a}/{prot_a_domain}", f"{prot_b}/{prot_b_domain}"))

    print("Skipped:", skips)
    domainG_graph = nx.Graph(domainG_random)
    print("Reassigned PPI/DDI graph:", domainG_graph)
    return domainG_graph
