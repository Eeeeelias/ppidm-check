import pickle
from networkx import expected_degree_graph
from networkx.classes.graphviews import generic_graph_view


def sort_dict(dictionary: dict):
    return {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1], reverse=True)}


gold_standard = pickle.load(open('gold_standard.pickle', 'rb'))

single_pfam = {}
for i, j in gold_standard:
    try:
        single_pfam[i] = single_pfam.get(i) + 1
    except TypeError:
        single_pfam[i] = 1

    try:
        single_pfam[j] = single_pfam.get(j) + 1
    except TypeError:
        single_pfam[j] = 1

single_pfam = sort_dict(single_pfam)

count = 0
for i,j in single_pfam.items():
    print(i, j)
    count += 1
    if count > 5:
        break


random_graph = expected_degree_graph(list(single_pfam.values()))

print(random_graph)
print(len(random_graph))
print(random_graph.nodes[1])