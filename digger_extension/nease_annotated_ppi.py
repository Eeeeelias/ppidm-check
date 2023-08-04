import pickle

import networkx as nx
import networkx.exception
import pandas as pd

path = "D:/programming/bachelor_projects/NEASE/nease/data"

ppi_graph: nx.Graph = pickle.load(open(path + "/network/PPI_Human.pkl", 'rb'))
print("PPIs:", ppi_graph)
ppis = [(i, j) if i < j else (j, i) for i, j in ppi_graph.edges]


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
    pdb_db: pd.DataFrame = pickle.load(open(pdb_path, 'rb'))
    print(pdb_db.columns)
    pdb_interactions = [(str(row['Interator gene 1']), str(row['Interator gene 2']))
                        if str(row['Interator gene 1']) < str(row['Interator gene 2'])
                        else (str(row['Interator gene 2']), str(row['Interator gene 1'])) for _, row in pdb_db.iterrows()]


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
    ppis_filtered = filter_by_ddi(path + "/network/graph_human_ext.pkl")
    ppis_filtered += filter_by_elm(path + "/database/ELM_interactions")
    # filter_by_pdb(path + "\\database\\pdb")
    pathways = pickle.load(open(path + "/pathways/pathways_human", 'rb'))
    filtered_graph = nx.Graph(ppis_filtered)
    print(filtered_graph)
    new_degrees = pathway_degree(pathways, filtered_graph)
    # pathways['Degree in the PPI/DDI_old'] = pathways['Degree in the PPI/DDI']
    # pathways['Degree in the structural PPI_old'] = pathways['Degree in the structural PPI']
    pathways['Degree in the structural PPI'] = new_degrees
    pickle.dump(pathways, open(path + "/pathways/pathways_human", 'wb'))

