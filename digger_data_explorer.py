import pickle
import networkx as nx
import pandas as pd
import mygene

from did_comparison import read_interactions, domine, did_2022, did_2017, inter_predicted

gid2name = pickle.load(open('data_digger/gid2name.pkl', 'rb'))


def load_graph(file_path: str):
    with open(file_path, 'rb') as g:
        ddi_G = pickle.load(g)
        print(ddi_G)
    ddi_set = set((e[0], e[1]) if e[0] < e[1] else (e[1], e[0]) for e in ddi_G.edges)
    return ddi_set


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


# make nice table
ddi_set = load_graph('data_digger/DomainG.pkl')
interaction = ('PF07525', 'PF00023')
domains_int = set()
for d1, d2 in ddi_set:
    domain1 = d1.split("/")[1]
    domain2 = d2.split("/")[1]
    if interaction == (domain1, domain2) or interaction == (domain2, domain1):
        print(d1, d2)
    domains_int.add((domain1, domain2))
print("Distinct domains:", len(domains_int))
# print(entrez_to_name_online('6874'))


ddi_pkl = read_interactions('data_digger/DDI.pkl_interactions.tsv')
domain_g = read_interactions('data_digger/DomainG.pkl_interactions.tsv')

print("DDI.pkl interactions:", len(set(ddi_pkl)))
print("DomainG.pkl interactions (DDI only, \"duplicates\" removed):", len(set(domain_g)))
domain_g_x = set(domain_g) - set(ddi_pkl)
print("Exclusively in DomainG.pkl:", len(domain_g_x))
# print(domain_g_x)
for name, ddi_source in zip(['domine', 'ddi_pkl', 'did_2017', 'did_2022', 'predicted'],
                            [domine, ddi_pkl, did_2017, did_2022, inter_predicted]):
    if ('PF07525', 'PF00023') in ddi_source:
        print(name, True)
