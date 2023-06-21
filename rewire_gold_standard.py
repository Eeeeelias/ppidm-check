import pickle

import networkx as nx
from networkx import expected_degree_graph
from networkx.classes.graphviews import generic_graph_view


def sort_dict(dictionary: dict):
    return {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1], reverse=True)}


def sorted_node_degrees(network: list[tuple]) -> dict:
    single_pfam = {}
    for i, j in network:
        try:
            single_pfam[i] = single_pfam.get(i) + 1
        except TypeError:
            single_pfam[i] = 1

        try:
            single_pfam[j] = single_pfam.get(j) + 1
        except TypeError:
            single_pfam[j] = 1

    return sort_dict(single_pfam)


gold_standard = pickle.load(open('train_set.pickle', 'rb'))

single_pfam = sorted_node_degrees(gold_standard)

pfam_ids = list(single_pfam.keys())
pfam_degrees = list(single_pfam.values())

# printing the domains with the highest edge degree
for i in range(5):
    print(pfam_ids[i], single_pfam[pfam_ids[i]])
print("")

# #######################
# making the random graph
# #######################

random_graph = expected_degree_graph(pfam_degrees, seed=42)

random_gold = set()
for i, j in random_graph.edges:
    random_gold.add((pfam_ids[i], pfam_ids[j]))

# printing the node degrees of the new graph
random_sorted = sorted_node_degrees(random_gold)
random_graph_ids = list(random_sorted.keys())

for i in range(5):
    print(random_graph_ids[i], random_sorted[random_graph_ids[i]])

# saving the graph for it to be used in PPIDM
pickle.dump(random_gold, open('random_train.pickle', 'wb'))

for key in single_pfam.keys():
    if key not in random_sorted.keys():
        random_sorted[key] = 0


# Visualising
import matplotlib.pyplot as plt

x = list(random_sorted.values())
y = list(single_pfam.values())
print(len(x), len(y))
plt.scatter(x, y)
# Adding labels to the plot
plt.xlabel('expected degree')
plt.ylabel('true degree')
plt.plot([0, x[0]], [0, y[0]], 'r--')
plt.title('Retaining node degree of the random network')

# Displaying the plot
plt.show()
plt.savefig("expected_node_degree.png")