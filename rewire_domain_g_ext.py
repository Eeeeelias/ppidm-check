import pickle
import networkx as nx
from matplotlib import pyplot as plt
from rewire_gold_standard import sorted_node_degrees

ext_graph: nx.Graph = pickle.load(open('data_digger/DomainG_extended.pkl', 'rb'))
print(ext_graph)
sorted_ext_graph = sorted_node_degrees(list(ext_graph.edges))
nodes = list(sorted_ext_graph.keys())
degrees = list(sorted_ext_graph.values())

random_ext_graph = nx.expected_degree_graph(degrees, seed=42)

node_mapping = {i: j for i, j in zip(range(len(nodes)), nodes)}
random_ext_graph: nx.Graph = nx.relabel_nodes(random_ext_graph, node_mapping)
print(random_ext_graph)
pickle.dump(random_ext_graph, open('data_nease/random_graph_human.pkl', 'wb'))

random_degrees = {}
random_interactions = list(random_ext_graph.edges)
for key in sorted_ext_graph.keys():
    random_degrees[key] = len(list(random_ext_graph.neighbors(key)))

# Visualising
plt.style.use('ggplot')
x = degrees
y = list(random_degrees.values())

plt.scatter(x, y, color='#1F77B4')
# Adding labels to the plot
plt.xlabel('expected degree')
plt.ylabel('true degree')
plt.plot([0, x[0]], [0, x[0]], 'r--')
plt.title('Node degree of the random extended DomainG network')

# Displaying the plot
# plt.savefig("pictures/domaing_expected_node_degree.png")
plt.show()
