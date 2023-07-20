import pickle
import networkx as nx
import pandas as pd
import mygene

from did_comparison import read_interactions, domine, did_2022, did_2017, inter_predicted

gid2name = pickle.load(open('data_digger/gid2name.pkl', 'rb'))


def load_graph(file_path: str, ppi_ddi=False):
    with open(file_path, 'rb') as g:
        ddi_G = pickle.load(g)
        print(ddi_G)
    if ppi_ddi:
        return set((e[0], e[1]) if e[0].split("/")[1] < e[1].split("/")[1] else (e[1], e[0]) for e in ddi_G.edges)
    return set((e[0], e[1]) if e[0] < e[1] else (e[1], e[0]) for e in ddi_G.edges)


def entrez_to_name_online(entrezID):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(entrezID, scopes='entrezgene', fields='symbol', species='human')
    return out


def export_table(file_path):
    # export nice table
    with open(file_path + '.tsv', 'w') as f:
        # header
        f.write('domain_1\tdomain_2\n')
        for ddi in domains_int:
            f.write(f'{ddi[0]}\t{ddi[1]}\n')


if __name__ == '__main__':
    # load ddi/ppi graph and add new domains to DomainG.pkl
    ddi_ppi_set = load_graph('data_digger/DomainG.pkl', ppi_ddi=True)
    domains_int = set()
    for d1, d2 in ddi_ppi_set:
        domain1 = d1.split("/")[1]
        domain2 = d2.split("/")[1]
        domains_int.add((domain1, domain2))
    print("Distinct domains:", len(domains_int))
    # print(entrez_to_name_online('6874'))
    # export_table("data_digger/DomainG.pkl_interactions")

    ddi_predicted = read_interactions('predicted_ddi_ppi.tsv.tsv', third_col='all')
    print("Predicted PPI/DDI interactions:", len(set(ddi_predicted)))
    for ddi in ddi_predicted:
        ddi_ppi_set.add(ddi)

    print("Exclusively in DomainG.pkl:", len(set(ddi_ppi_set) - set(ddi_predicted)))

    ext_graph = nx.Graph(ddi_ppi_set)
    print(ext_graph)
    # pickle.dump(ext_graph, open("data_digger/DomainG_extended.pkl", 'wb'))
    # Gold: 384.737
    # Silver: 109.185
    # Bronze: 565.520

    # extending the DDI graph
    ddi_set = load_graph('data_digger/DDI.pkl')
    print(list(set(inter_predicted) - set(ddi_set))[:5])
    for interaction in inter_predicted:
        ddi_set.add(interaction)
    ext_ddi_graph = nx.Graph(ddi_set)
    print(ext_ddi_graph)
    # pickle.dump(ext_ddi_graph, open("data_digger/DDI_extended.pkl", 'wb'))
    interaction = ('PF00640', 'PF07714')
