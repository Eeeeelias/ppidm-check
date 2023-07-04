import pickle

import networkx as nx

g = nx.read_edgelist('digger_ddis_annotated.tsv', delimiter='\t', comments='d',
                     data=[('anno_by_3did', int), ('anno_by_domine', int)])

print(g)

with open('resultdata/result-all', 'r') as f:
    for line in f.readlines():
        if not line.startswith("PF"):
            continue
        parts = line.split("\t")
        g.add_edge(parts[0], parts[1])

print(g)

domain_g = pickle.load(open('data_digger/DomainG.pkl', 'rb'))
print(domain_g)

combination: nx.Graph = nx.compose(g, domain_g)
print(combination)

ddi = pickle.load(open('data_digger/DDI.pkl', 'rb'))
print(ddi)
