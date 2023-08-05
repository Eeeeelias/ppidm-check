import pickle
import networkx as nx
import pandas as pd

path = "D:/programming/bachelor_projects/NEASE/nease/data"

ppi_graph: nx.Graph = pickle.load(open(path + "/network/PPI_Human.pkl", 'rb'))
print("PPIs:", ppi_graph)
ppis = [(i, j) if i < j else (j, i) for i, j in ppi_graph.edges]


def extract_pdb_supported_ppis():
    pdb_info = pd.read_csv("D:/Downloads/PPI_table_proteome_95_seq_id.tsv", sep="\t", low_memory=False)
    mart_info = pd.read_csv("sourcedata/mart_export.txt", sep="\t")

    uniprot_entrez = dict()
    for idx, row in mart_info.iterrows():
        if pd.isna(row['NCBI gene (formerly Entrezgene) ID']):
            continue
        uniprot_entrez[row['UniProtKB/Swiss-Prot ID']] = row['NCBI gene (formerly Entrezgene) ID']

    print("Proteins mapped to entrez IDs:", len(uniprot_entrez))

    interactions = set()
    pdb_len = len(pdb_info)
    for idx, row in pdb_info.iterrows():
        # just for progress
        if idx % 100_000 == 0:
            print(round(idx / pdb_len, 2) * 100, "%")
        try:
            # get entrez id of uniprot protein
            gene1 = str(int(uniprot_entrez.get(row['u_ac_1'])))
            gene2 = str(int(uniprot_entrez.get(row['u_ac_2'])))
        except KeyError:
            continue
        except TypeError:
            continue
        # add it in sorted order, so it's easier to compare later
        interactions.add((gene1, gene2)) if gene1 < gene2 else interactions.add((gene2, gene1))

    # write it to a file
    print("PDB interactions:", len(interactions))
    with open("residue_supported_interactions.tsv", 'w') as f:
        for i, j in interactions:
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
    ppis_filtered = filter_by_ddi(path + "/network/graph_human_old.pkl")
    ppis_filtered += filter_by_elm(path + "/database/ELM_interactions")
    ppis_filtered += filter_by_pdb("../residue_supported_interactions_exons.tsv")
    pathways = pickle.load(open(path + "/pathways/pathways_human", 'rb'))
    filtered_graph = nx.Graph(ppis_filtered)
    print(filtered_graph)
    # new_degrees = pathway_degree(pathways, filtered_graph)
    # # pathways['Degree in the PPI/DDI_old'] = pathways['Degree in the PPI/DDI']
    # # pathways['Degree in the structural PPI_old'] = pathways['Degree in the structural PPI']
    # pathways['Degree in the structural PPI'] = new_degrees
    # pickle.dump(pathways, open(path + "/pathways/pathways_human", 'wb'))

