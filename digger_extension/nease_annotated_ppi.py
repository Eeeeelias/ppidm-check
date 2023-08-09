import pickle
import networkx as nx
import pandas as pd

path = "D:/programming/bachelor_projects/NEASE/nease/data"

ppi_graph: nx.Graph = pickle.load(open(path + "/network/PPI_Human.pkl", 'rb'))
print("PPIs:", ppi_graph)
ppis = [(i, j) if i < j else (j, i) for i, j in ppi_graph.edges]


def extract_pdb_supported_ppis():
    pdb_info = pd.read_csv("../data_nease/PPI_table_clean.tsv", sep="\t")
    df = pd.read_csv("../data_nease/HUMAN_9606_idmapping_GeneIDs.dat", sep="\t", header=None)
    df.columns = ['UniProtID', "GeneID", "EntrezID"]
    id_mapping = df.set_index('UniProtID')['EntrezID'].to_dict()
    # mart_info = pd.read_csv("../sourcedata/mart_export.txt", sep="\t")
    # id_mapping = mart_info.set_index('UniProtKB/Swiss-Prot ID')['NCBI gene (formerly Entrezgene) ID'].to_dict()

    print("Proteins mapped to entrez IDs:", len(id_mapping))
    orig_interactions = set()
    interactions = set()
    key_errors = 0
    for idx, row in pdb_info.iterrows():
        try:
            orig_interactions.add((row['u_ac_1'], row['u_ac_2']))
            # get entrez id of uniprot protein
            gene1 = str(id_mapping.get(row['u_ac_1']))
            gene2 = str(id_mapping.get(row['u_ac_2']))
        except KeyError:
            key_errors += 1
            continue
        if gene1 == 'None' or gene2 == 'None':
            continue
        # add it in sorted order, so it's easier to compare later
        interactions.add((gene1, gene2)) # if gene1 < gene2 else interactions.add((gene2, gene1))

    print("Key errors:", key_errors)
    # write it to a file
    print("PDB interactions:", len(interactions))
    with open("../residue_supported_interactions.tsv", 'w') as f:
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
    pdb_db_reversed = pdb_db
    for i, j in pdb_db:
        pdb_db_reversed.add((j, i))
    filtered = set(ppis) & pdb_db_reversed
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


def combine_exon_residue_interactions(exon_level: str, residue_level: str):
    exon_interactions = set()
    residue_interactions = set()

    with open(exon_level, 'r') as f:
        for line in f.readlines():
            el1 = line.split("\t")[0]
            el2 = line.split("\t")[1].strip()
            exon_interactions.add((el1, el2))
    with open(residue_level, 'r') as f:
        for line in f.readlines():
            el1 = line.split("\t")[0]
            el2 = line.split("\t")[1].strip()
            residue_interactions.add((el1, el2))

    print("Exon level interactions:", len(exon_interactions))
    print("Residue level interactions:", len(residue_interactions))
    print("Combined:", len(exon_interactions | residue_interactions))


if __name__ == '__main__':
    extract_pdb_supported_ppis()
    ppis_filtered = filter_by_ddi(path + "/network/graph_human_old.pkl")
    ppis_filtered += filter_by_elm(path + "/database/ELM_interactions")
    ppis_filtered += filter_by_pdb("../residue_supported_interactions.tsv")
    pathways = pickle.load(open(path + "/pathways/pathways_human", 'rb'))
    filtered_graph = nx.Graph(ppis_filtered)
    print("PPI graph supported by DDI/ELM/PDB:", filtered_graph)
    # new_degrees = pathway_degree(pathways, filtered_graph)
    # # pathways['Degree in the PPI/DDI_old'] = pathways['Degree in the PPI/DDI']
    # # pathways['Degree in the structural PPI_old'] = pathways['Degree in the structural PPI']
    # pathways['Degree in the structural PPI'] = new_degrees
    # pickle.dump(pathways, open(path + "/pathways/pathways_human", 'wb'))

