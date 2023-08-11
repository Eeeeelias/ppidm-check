import pickle
import networkx as nx
import numpy as np
import pandas as pd

path = "D:/programming/bachelor_projects/NEASE/nease/data"

ppi_graph: nx.Graph = pickle.load(open(path + "/network/PPI_Human.pkl", 'rb'))
print("PPIs:", ppi_graph)
ppis = [(i, j) if i < j else (j, i) for i, j in ppi_graph.edges]


def combine_uniprot_mappings():
    # Read UniProt to Entrez mapping
    gene_id_mapping = pd.read_csv("../data_nease/HUMAN_9606_idmapping_GeneIDs.dat", sep="\t", header=None)
    gene_id_mapping.columns = ['UniProtID', 'GeneID', 'EntrezID']
    gene_id_mapping['EntrezID'] = gene_id_mapping['EntrezID'].astype(str)
    id_mapping = gene_id_mapping.set_index('UniProtID')['EntrezID'].to_dict()

    # Read UniProt to Ensembl mapping
    ensembl_info = pd.read_csv("../data_nease/uniprot_ensembl_mapping.tsv", sep="\t", header=None)
    ensembl_info.columns = ['UniProtID', 'Ensembl', 'EnsemblID']
    ensembl_mapping = ensembl_info.set_index('UniProtID')['EnsemblID'].to_dict()

    # Read Ensembl to Entrez mapping
    mart_info = pd.read_csv("../sourcedata/mart_export.txt", sep="\t",
                            dtype={'Gene stable ID': str, 'UniProtKB/Swiss-Prot ID': str,
                                   'NCBI gene (formerly Entrezgene) ID': str})
    mart_mapping = mart_info.set_index('Gene stable ID')['NCBI gene (formerly Entrezgene) ID'].to_dict()
    mart_mapping_entrez = mart_info.set_index('UniProtKB/Swiss-Prot ID')['NCBI gene (formerly Entrezgene) ID'].to_dict()
    mart_mapping_entrez.update(id_mapping)
    id_mapping = mart_mapping_entrez

    for i, j in ensembl_mapping.items():
        # skip most of this
        if id_mapping.get(i) is not None:
            continue

        ens_map = mart_mapping.get(j)
        # skip if no ensembl map can be found
        if ens_map is None:
            continue
        id_mapping[i] = ens_map

    tmp_mapping = dict()
    for i, j in id_mapping.items():
        if type(j) != str:
            continue
        tmp_mapping[i] = j

    return tmp_mapping


def extract_pdb_supported_ppis():
    # Read protein interactions data
    # download the id mapping and extract the necessary geneIds like so:
    # cat PPI_table_proteome_95_seq_id.tsv | cut -d$'\t' -f1,6 | sed 's/-[0-9]//g' | sort -r | uniq >
    # PPI_table_clean.tsv
    pdb_info = pd.read_csv("../data_nease/PPI_table_clean.tsv", sep="\t")
    id_mapping = combine_uniprot_mappings()

    interactions = set()
    interactions_n_r = set()
    nones = 0

    for _, row in pdb_info.iterrows():
        try:
            gene1 = id_mapping.get(row['u_ac_1'])
            gene2 = id_mapping.get(row['u_ac_2'])
        except KeyError:
            continue

        if gene1 is None or gene2 is None:
            nones += 1
            continue

        interactions.add((gene1, gene2))
        gene1, gene2 = sorted([gene1, gene2])
        interactions_n_r.add((gene1, gene2))

    print("Proteins mapped to Entrez IDs:", len(id_mapping))
    print("Interactions not mapped:", nones)
    print("PDB interactions:", len(interactions))
    print("PDB non-redundant interactions:", len(interactions_n_r))

    # Write interactions to a file
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

