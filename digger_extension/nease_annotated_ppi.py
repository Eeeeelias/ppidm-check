import pickle

import mygene
import networkx as nx
import pandas as pd
import rewire_domain_g_ext as rg

path = "D:/programming/bachelor_projects/NEASE/nease/data"

ppi_graph: nx.Graph = pickle.load(open(path + "/network/PPI_Human.pkl", 'rb'))
print("PPIs:", ppi_graph)
ppis = [(i, j) if i < j else (j, i) for i, j in ppi_graph.edges]


def name_to_entrez(gene, mapping):
    try:
        name = mapping[mapping['Gene name'] == str(gene)]['NCBI gene ID'].unique()

        if len(name) == 0:
            # try to convert it online
            mg = mygene.MyGeneInfo()
            out = mg.query(f"symbol:{gene}", fields='entrezgene', species="human")
            return out['hits'][0]['_id']

        return name[0]
    except:
        return None


def extract_pdb_interactions():
    pdb = pickle.load(open("D:\programming\\bachelor_projects\\NEASE\\nease\data\database\pdb", 'rb'))
    human = pickle.load(open("D:\\programming\\bachelor_projects\\NEASE\\nease\\data\\database\\Human", 'rb'))

    # mapping_dict = {gene: name_to_entrez(gene, human) for gene in pdb['symbol'].unique()}
    # pickle.dump(mapping_dict, open("../tmp_mapping_dict.pkl", 'wb'))

    mapping_dict = pickle.load(open("../tmp_mapping_dict.pkl", 'rb'))

    interactions = set()
    for _, row in pdb.iterrows():
        entrez_retrieved = mapping_dict.get(row['symbol'])
        if entrez_retrieved is None or str(row['entrezgene']) == 'nan':
            continue
        interactions.add((str(entrez_retrieved), str(row['entrezgene'])))

    print(len(interactions))
    G = nx.Graph(interactions)
    print(G)

    interaction_set = set()
    for i, j in interactions:
        interaction_set.add((i, j)) if i < j else interaction_set.add((j, i))

    with open("residue_supported_interactions.tsv", 'w') as f:
        for i, j in interaction_set:
            f.write(f"{i}\t{j}\n")


def filter_by_ddi(ddi_path):
    ddi_graph: nx.Graph = pickle.load(open(ddi_path, 'rb'))
    print("DDIs:", ddi_graph)

    ddis = [(i.split("/")[0], j.split("/")[0]) if i.split("/")[0] < j.split("/")[0] else
            (j.split("/")[0], i.split("/")[0]) for i, j in ddi_graph.edges]


    filtered = set(ddis) & set(ppis)
    print(f"{len(filtered):,} interactions supported by DDI")
    return list(filtered)


def filter_by_elm(elm_path):
    elm_df: pd.DataFrame = pickle.load(open(elm_path, 'rb'))
    elm_interators = [(str(row['Interator gene 1']), str(row['Interator gene 2']))
                      if str(row['Interator gene 1']) < str(row['Interator gene 2'])
                      else (str(row['Interator gene 2']), str(row['Interator gene 1'])) for _, row in elm_df.iterrows()]

    filtered = set(elm_interators) & set(ppis)
    print(f"{len(filtered):,} interactions supported by ELM")
    return list(filtered)


def filter_by_pdb(pdb_path):
    pdb_db = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            split = line.strip().split("\t")
            pdb_db.add((split[0], split[1]))
    filtered = set(ppis) & pdb_db
    print(f"{len(filtered):,} interactions supported by PDB")
    return list(filtered)


def pathway_degree(pathways: pd.DataFrame, ppi_ddi_graph: nx.Graph):
    new_degrees = []
    for value in pathways['entrez_gene_ids']:
        current_pathway = 0
        for gene in value:
            try:
                current_pathway += ppi_ddi_graph.degree[str(gene)]
            except KeyError:
                continue
        new_degrees.append(current_pathway)
    return new_degrees


if __name__ == '__main__':
    extended_graph = "C:\\Users\\gooog\\miniconda3\\Lib\\site-packages\\nease\\data\\network\\graph_human_ext.pkl"
    # extract_pdb_interactions()
    ppis_filtered = filter_by_ddi(extended_graph)
    ppis_filtered += filter_by_elm(path + "/database/ELM_interactions")
    ppis_filtered += filter_by_pdb("../residue_supported_interactions.tsv")
    filtered_graph = nx.Graph(ppis_filtered)
    print("PPI graph supported by DDI/ELM/PDB:", filtered_graph)

    # this is only relevant for the random graph
    random_supported_graph, degrees, random_degrees = rg.rewire_graph(filtered_graph)
    rg.visualize(degrees, random_degrees, title="Retaining node degree in the random PPI graph",
                 file_path="../pictures/domaing_expected_node_degree.png")
    random_ppi_ddi_graph = rg.reassign_ddis(pickle.load(open(extended_graph, 'rb')), random_supported_graph)
    pickle.dump(random_ppi_ddi_graph, open("C:\\Users\\gooog\\miniconda3\\Lib\\site-packages\\nease\\data\\network\\"
                                     "graph_human_ext_random.pkl", 'wb'))

    # continue here for normal execution
    pathways = pickle.load(open(path + "/pathways/pathways_human", 'rb'))
    new_degrees = pathway_degree(pathways, random_supported_graph)
    # pathways['Degree in the PPI/DDI_old'] = pathways['Degree in the PPI/DDI']
    # pathways['Degree in the structural PPI_old'] = pathways['Degree in the structural PPI']
    pathways['Degree in the structural PPI'] = new_degrees
    pickle.dump(pathways, open(path + "/pathways/pathways_human_ext_random", 'wb'))
