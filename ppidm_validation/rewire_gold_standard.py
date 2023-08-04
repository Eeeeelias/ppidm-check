import pickle
from networkx import expected_degree_graph
import matplotlib.pyplot as plt


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


def change_random_ordering(real_order: list, random_degrees: dict) -> dict:
    return {domain: random_degrees[domain] for domain in real_order}


if __name__ == '__main__':
    gold_standard = pickle.load(open('../pickles/train_set.pickle', 'rb'))

    single_pfam = sorted_node_degrees(gold_standard)

    pfam_ids = list(single_pfam.keys())
    pfam_degrees = list(single_pfam.values())


    # #######################
    # making the random graph
    # #######################

    random_graph = expected_degree_graph(pfam_degrees, seed=42)

    random_gold = set()
    for i, j in random_graph.edges:
        random_gold.add((pfam_ids[i], pfam_ids[j]))

    # printing the node degrees of the new graph
    random_sorted = sorted_node_degrees(random_gold)

    # saving the graph for it to be used in PPIDM
    pickle.dump(random_gold, open('../pickles/random_train.pickle', 'wb'))

    # adding zeroes just for visualisation
    for key in single_pfam.keys():
        if key not in random_sorted.keys():
            random_sorted[key] = 0

    random_sorted = change_random_ordering(pfam_ids, random_sorted)
    random_graph_ids = list(random_sorted.keys())

    # printing the domains with the highest edge degree
    for i in range(5):
        print(pfam_ids[i], single_pfam[pfam_ids[i]])
    print("")

    for i in range(5):
        print(random_graph_ids[i], random_sorted[random_graph_ids[i]])

    # Visualising
    plt.style.use('ggplot')
    x = list(random_sorted.values())
    y = list(single_pfam.values())
    plt.scatter(x, y, color='#1F77B4')
    # Adding labels to the plot
    plt.xlabel('expected degree')
    plt.ylabel('true degree')
    plt.plot([0, x[0]], [0, x[0]], 'r--')
    plt.title('Retaining node degree of the random network')

    # Displaying the plot
    plt.savefig("../pictures/expected_node_degree.png")
    plt.show()
