import pickle

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2


def read_interactions(file: str):
    interactions = []
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.strip().split("\t")
            interactions.append((line[0], line[1]))
    return interactions


def venn_diagrams(did_2017, did_2022, predicted, domine, category='gold'):
    fig, (ax1, ax2) = plt.subplots(1, 2)

    venn3([set(did_2017), set(predicted), set(did_2022)], set_labels=('3did2017', 'predicted', '3did2022',), ax=ax1)
    venn2([set(domine), set(predicted)], set_labels=('domine', 'predicted'), ax=ax2)

    plt.suptitle(f'Overlap known new databases and predicted ({category}) interactions')
    plt.savefig(f'overlap_domine_3did_{category}.png')
    plt.show()


def remove_unknown_domains(known_domains, new_domains):
    clean = []
    for i, j in new_domains:
        if i in known_domains and j in known_domains:
            clean.append((i, j))
    return clean


# loading all files
did_2022 = read_interactions('resultdata/3did_2022')
did_2017 = read_interactions('resultdata/3did')
domine = read_interactions('resultdata/domine')
interactions_gold = read_interactions('resultdata/interactions_gold')
inter_silver = read_interactions("resultdata/interactions_silver")
inter_bronze = read_interactions("resultdata/interactions_bronze")
inter_predicted = interactions_gold + inter_silver + inter_bronze
all_known_ids = pickle.load(open('all_known_ids.pickle', 'rb'))

# removing domains that were not known at the time of PPIDM
did_2022_clean = remove_unknown_domains(all_known_ids, did_2022)
domine_clean = remove_unknown_domains(all_known_ids, domine)


# visualisation
venn_diagrams(did_2017, did_2022_clean, inter_predicted, domine_clean, category='all')

print("Length of interactions_gold:", len(inter_predicted))
print("Length of did_2022:", len(did_2022))
print("Length of did_2017:", len(did_2017))
print("Unique 2022:", len(set(did_2017) - set(did_2022)))
print("Interaction overlap gold and did_2017:", len(set(interactions_gold) & set(did_2017)))
print("Interaction overlap gold and did_2022:", len(set(interactions_gold) & set(did_2022)))
