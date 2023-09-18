import pickle

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3, venn2
from upsetplot import from_contents, UpSet
import seaborn as sns
from ppidm_run.main import source_address

plt.style.use('ggplot')


def sort_dict(dictionary: dict):
    if isinstance(list(dictionary.values())[0], list):
        return {k: v for k, v in sorted(dictionary.items(), key=lambda item: sum(item[1]), reverse=True)}
    return {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1], reverse=True)}


def read_interactions(file: str, third_col=None):
    interactions = []
    header = True
    with open(file, 'r') as f:
        for line in f.readlines():
            if not header:
                header = True
                continue
            line = line.strip().split("\t")
            if third_col is not None:
                if line[2] not in third_col:
                    continue
            assoc = (line[0], line[1]) if line[0] < line[1] else (line[1], line[0])
            # assoc = (line[0], line[1])
            interactions.append(assoc)
    return interactions


def filter_domine():
    with open(source_address + 'INTERACTION.txt', 'r') as f:
        associations = []
        for line in f.readlines():
            if 'PF' not in line:
                continue
            line_sp = line.split('|')
            if line_sp[17] == 'HC' or line_sp[17] == 'MC':
                PF1 = line_sp[0]
                PF2 = line_sp[1]
                if PF1 < PF2:
                    associations.append((PF1, PF2))
                else:
                    associations.append((PF2, PF1))
    return associations


def filter_domine_sources():
    sources = {2: 'iPfam', 3: '3did', 4: 'ME', 5: 'RCDP', 6: 'P-value', 7: 'Interdom', 8: 'DPEA', 9: 'PE',
               10: 'GPE', 11: 'DIPD', 12: 'RDFF', 13: 'K-GIDDI', 14: 'Insite', 15: 'DomainGA', 16: 'DIMA'}
    interactions = {src: [] for src in list(sources.values())}
    with open(source_address + 'INTERACTION.txt', 'r') as f:
        for line in f.readlines():
            if 'PF' not in line:
                continue
            line_sp = line.split('|')
            PF1 = ""
            PF2 = ""
            for index, part in enumerate(line_sp):
                if index >= 17:
                    break
                if index == 0:
                    PF1 = part
                    continue
                if index == 1:
                    PF2 = part
                    continue
                # assoc = f"{PF1};{PF2}" if PF1 > PF2 else f"{PF2};{PF1}"
                assoc = (PF1, PF2) if PF1 < PF2 else (PF2, PF1)
                if part == "1":
                    interactions[sources[index]].append(assoc)
    return interactions


def remove_unknown_domains(known_domains, new_domains):
    clean = []
    for i, j in new_domains:
        if i in known_domains and j in known_domains:
            clean.append((i, j))
    return clean


def venn_diagrams(did_2017, did_2022, predicted, domine, category='gold'):
    fig, (ax1, ax2) = plt.subplots(1, 2)

    venn3([set(did_2017), set(predicted), set(did_2022)], set_labels=('3did2017', 'predicted', '3did2022',), ax=ax1)
    venn2([set(domine), set(predicted)], set_labels=('domine', 'predicted'), ax=ax2)

    plt.suptitle(f'Overlap known new databases and predicted ({category}) interactions')
    plt.savefig(f'../pictures/overlap_domine_3did_{category}.png')
    plt.show()


def upset_plots(unformatted_dict: dict, title='Generic UpsetPlot', min_subset_size=None):
    upset_format = from_contents(unformatted_dict)
    if min_subset_size:
        ax_dict = UpSet(upset_format, show_counts=True, min_subset_size=min_subset_size).plot()
    else:
        ax_dict = UpSet(upset_format, show_counts=True).plot()
    plt.savefig(f'../pictures/{title}.png')
    plt.show()


def domine_pairwise_comparison(predicted: list[tuple[str, str]], domine: dict[list[tuple[str, str]]], save=True):
    bold = '\033[1m'
    end = '\033[0m'
    overlap_relative_domine = {}
    for key, value in domine.items():
        overlap = len(set(value) & set(predicted))
        relative_overlap = overlap / len(set(value))
        overlap_relative_domine[key] = relative_overlap * 100
        print(f"Predicted and {bold}{key}{end} overlap for {overlap}/{len(set(value))} "
              f"({round(relative_overlap * 100, 1)}%) interactions")

    overlap_relative_domine = sort_dict(overlap_relative_domine)
    x = list(overlap_relative_domine.keys())
    y = list(overlap_relative_domine.values())
    plt.figure(figsize=(8, 6))

    plt.bar(x, y, color='black')
    plt.xlabel('Sources')
    plt.ylabel('Overlap (in %)')
    plt.title('Overlap DOMINE sources with predicted interactions')

    for i, v in enumerate(y):
        plt.text(i, v, str(round(v, 1)), ha='center', va='bottom')

    plt.xticks(rotation=45)
    plt.tight_layout()
    if save:
        plt.savefig('../pictures/domine_comparison.png')
    plt.show()


def domine_pairwise_comparison_categories(predicted: dict[str, list[tuple[str, str]]],
                                              domine: dict[str, list[tuple[str, str]]], save=True):
    overlap_relative_domine = {}
    for key, value in domine.items():
        relative_overlaps = []
        for pred_key in predicted.keys():
            overlap = len(set(value) & set(predicted[pred_key]))
            relative_overlaps.append((overlap / len(set(value))) * 100)
        overlap_relative_domine[key] = relative_overlaps
        print(f"Source: {key} \t Overlaps: {relative_overlaps} \t Sum: {sum(relative_overlaps)}")

    overlap_relative_domine = sort_dict(overlap_relative_domine)
    df = pd.DataFrame(overlap_relative_domine)

    df.columns = df.columns.str.replace('-', '')
    df = df.add_prefix('Overlap_')
    df.loc[3] = df.sum()
    df['category'] = ['gold', 'silver', 'bronze', 'total']
    df = pd.wide_to_long(df, ['Overlap'], i='category', j='source', sep='_', suffix='(\d+|\w+)')
    df = df.reset_index()
    sns.set(rc={'figure.figsize': (14, 4)})
    sns.barplot(data=df, x='source', y='Overlap', hue='category', hue_order=['total', 'gold', 'silver', 'bronze'])

    # plt.tight_layout()
    # plt.xticks(rotation=45)
    plt.xlabel("Sources")
    plt.ylabel("Overlap (in %)")
    plt.title("Overlap DOMINE sources with predicted by category")

    if save:
        plt.savefig('../pictures/domine_comparison_category.png')
    plt.show()


# loading all files
did_2022 = read_interactions('../resultdata/3did_2022')
did_2017 = read_interactions('../resultdata/3did')
domine = filter_domine()
# domine = filter_domine_sources()
interactions_gold = read_interactions('../resultdata/interactions_gold')
inter_silver = read_interactions("../resultdata/interactions_silver")
inter_bronze = read_interactions("../resultdata/interactions_bronze")
inter_predicted = interactions_gold + inter_silver + inter_bronze
inter_predicted_dict = {'gold': interactions_gold, 'silver': inter_silver, 'bronze': inter_bronze}
all_known_ids = pickle.load(open('../pickles/all_known_ids.pickle', 'rb'))


if __name__ == '__main__':
    # removing domains that were not known at the time of PPIDM
    did_2022_clean = remove_unknown_domains(all_known_ids, did_2022)

    # visualisation
    # domine_pairwise_comparison_categories(predicted=inter_predicted_dict, domine=filter_domine_sources(), save=True)

    # did = {'did_2017': did_2017, 'did_2022': did_2022_clean, 'predicted': inter_predicted}
    # domine['predicted'] = inter_predicted
    # upset_plots(did, 'upset_did_comparison_count')
    # upset_plots(domine, 'upset_domine_comparison', min_subset_size=100)
    # venn_diagrams(did_2017, did_2022_clean, inter_predicted, domine, category='all')


    # random info

    # print("overlap 3did_2022:", (1 - (len(set(did_2022) - set(did_2017) - set(inter_predicted)) /
    # len(set(did_2022) - set(did_2017)))) * 100, "%")
    # print(len(interactions_gold), len(inter_silver), len(inter_bronze))
    print("domine:", len(domine))
    print("total predicted:", len(set(domine) & set(inter_predicted)))
    print("Length of interactions_gold:", len(inter_predicted))
    # print("Length of did_2022:", len(did_2022))
    # print("Length of did_2017:", len(did_2017))
    # print("Unique 2022:", len(set(did_2017) - set(did_2022)))
    # print("Interaction overlap gold and did_2017:", len(set(interactions_gold) & set(did_2017)))
    # print("Interaction overlap gold and did_2022:", len(set(interactions_gold) & set(did_2022)))
