import glob
import multiprocessing
import pickle
import threading

from tqdm import tqdm

from did_comparison import inter_predicted
import networkx as nx
import pandas as pd


def read_protein_interactions(file_path: str, predicted_int: set[tuple], uniprot_to_entrez: pd.DataFrame):
    interactions = set()
    with open(file_path, 'r') as f:
        all_lines = f.readlines()
        for line in tqdm(all_lines, desc=f'{file_path[11:]} Progress'):
            l_split = line.strip().split("\t")
            domain1 = l_split[0]
            domain2 = l_split[2]
            if '_' in domain1 or '_' in domain2:
                continue
            if (domain1, domain2) not in predicted_int and (domain2, domain1) not in predicted_int:
                continue
            prot1 = l_split[1].split('_')[0]
            prot2 = l_split[1].split('_')[1]
            prot1 = uniprot_to_entrez.loc[uniprot_to_entrez['UniProtID'] == prot1, 'EntrezID']
            prot2 = uniprot_to_entrez.loc[uniprot_to_entrez['UniProtID'] == prot2, 'EntrezID']
            if prot1.empty or prot2.empty:
                continue
            else:
                prot1 = prot1.iloc[0]
                prot2 = prot2.iloc[0]
            interactions.add((f"{prot1}/{domain1}", f"{prot2}/{domain2}"))
    print(f"Read {len(interactions)} interactions in total")
    return interactions


def im_dying(file_path: str):
    return read_protein_interactions(file_path, inter_predicted, mart_table)

# g = nx.read_edgelist('digger_ddis_annotated.tsv', delimiter='\t', comments='d',
#                      data=[('anno_by_3did', int), ('anno_by_domine', int)])
#
# print(g)
#
# with open('resultdata/result-all', 'r') as f:
#     for line in f.readlines():
#         if not line.startswith("PF"):
#             continue
#         parts = line.split("\t")
#         g.add_edge(parts[0], parts[1])


mart_table = pd.read_csv('sourcedata/mart_export.txt', sep='\t', dtype={'Gene stable ID': str,
                                                                        'UniProtKB/Swiss-Prot ID': str,
                                                                        'NCBI gene (formerly Entrezgene) ID': str})
mart_table.rename(columns={'Gene stable ID': 'GeneID', 'UniProtKB/Swiss-Prot ID': 'UniProtID',
                           'NCBI gene (formerly Entrezgene) ID': 'EntrezID'}, inplace=True)


def parallel_execution(inputs):
    pool = multiprocessing.Pool()
    results = pool.map(im_dying, inputs)
    pool.close()
    pool.join()

    return results


if __name__ == '__main__':
    input_list = ["resultdata/source" + source + "pfam_ordered" for source in
                  ["1_intact", "2_mint", "3_dip", "4_biogrid", "6_sifts", "5_string-exp", "5_string-rest", ]]
    output_list = parallel_execution(input_list[0:3])
    output_list += parallel_execution(input_list[3:5])
    output_list += parallel_execution(input_list[5:7])

    all_interactions = set()
    for s in output_list:
        all_interactions = all_interactions.union(s)

    with open("predicted_ddi_ppi" + '.tsv', 'w') as f:
        for ddi in all_interactions:
            f.write(f'{ddi[0]}\t{ddi[1]}\n')
