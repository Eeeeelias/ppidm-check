import multiprocessing
import re
import timeit

from ppidm_validation.did_comparison import inter_predicted
import pandas as pd


def read_protein_interactions(file_path: str, predicted_int: set[tuple], uniprot_to_entrez: pd.DataFrame):
    interactions = set()
    uniprot_to_entrez_dict = uniprot_to_entrez.set_index('UniProtID')['EntrezID'].to_dict()
    start = timeit.default_timer()
    pattern = re.compile(r"(PF.*\d)\t(.*)_(.*)\t(PF.*)")
    counter = 0
    with open(file_path, 'r') as f:
        print(f"Reading: {file_path}")
        for line in f:
            counter += 1
            if counter % 100000 == 0:
                print(f"{file_path}: {counter:,} lines")
            match = pattern.match(line.strip())
            try:
                domain1, prot1, prot2, domain2 = match.groups()
            except AttributeError:
                continue

            if '_' in domain1 or '_' in domain2:
                continue
            if (domain1, domain2) not in predicted_int and (domain2, domain1) not in predicted_int:
                continue
            prot1 = uniprot_to_entrez_dict.get(prot1)
            prot2 = uniprot_to_entrez_dict.get(prot2)
            if prot1 is None or prot2 is None:
                continue
            interactions.add((f"{prot1}/{domain1}", f"{prot2}/{domain2}"))

    print(f"Read {len(interactions)} interactions in {round((timeit.default_timer() - start) / 60 / 60, 1)} hours from "
          f"{file_path}")
    return interactions


def _wrapper_read_protein_interactions(file_path: str):
    return read_protein_interactions(file_path, inter_predicted, mart_table)


def add_classification(file_path: str, classifications_info: str):
    print("Reading the current PPI/DDI interactions")
    interactions = []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            domain1_ppi = line.split("\t")[0]
            domain2_ppi = line.strip().split("\t")[1]
            domain1 = domain1_ppi.split("/")[1]
            domain2 = domain2_ppi.split("/")[1]
            acc = (domain1_ppi, domain2_ppi) if domain1 < domain2 else (domain2_ppi, domain1_ppi)
            interactions.append(acc)

    print("Reading all the DDI classes")
    output_ddis = []
    tmp = set()
    classes = dict()
    with open(classifications_info, 'r') as f:
        for line in f.readlines():
            if not line.startswith("PF"):
                continue
            line_split = line.strip().split("\t")
            domain1 = line_split[0]
            domain2 = line_split[1]
            ddi_class = line_split[19]
            classes[(domain1, domain2)] = ddi_class

    for i, j in interactions:
        output_ddis.append((i, j, classes[(i.split("/")[1], j.split("/")[1])]))
        tmp.add((i, j))

    if len(output_ddis) != len(interactions):
        print(f"Something has gone wrong, not all DDIs in output anymore ({len(output_ddis)}/{len(interactions)}), "
              f"aborting")
        print(list(set(interactions) - set(tmp))[:5])
        return
    print("Writing the classes to file")
    with open(file_path + ".tsv", 'w') as f:
        f.write('domain_1\tdomain_2\tclass\n')
        for interact in output_ddis:
            f.write(f"{interact[0]}\t{interact[1]}\t{interact[2]}\n")


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


mart_table = pd.read_csv('../sourcedata/mart_export.txt', sep='\t', dtype={'Gene stable ID': str,
                                                                        'UniProtKB/Swiss-Prot ID': str,
                                                                        'NCBI gene (formerly Entrezgene) ID': str})
mart_table.rename(columns={'Gene stable ID': 'GeneID', 'UniProtKB/Swiss-Prot ID': 'UniProtID',
                           'NCBI gene (formerly Entrezgene) ID': 'EntrezID'}, inplace=True)


def parallel_execution(inputs):
    pool = multiprocessing.Pool()
    results = pool.map(_wrapper_read_protein_interactions, inputs)
    pool.close()
    pool.join()

    return results


def write_ddi_ppi_connection():
    input_list = ["resultdata/source" + source + "pfam_ordered" for source in
                  ["7_hprd", "2_mint", "3_dip", "1_intact", "4_biogrid", "6_sifts", "5_string-exp", "5_string-rest", ]]
    output_list = parallel_execution(input_list)
    # output_list += parallel_execution(input_list[3:5])
    # output_list += parallel_execution(input_list[5:7])

    all_interactions = set()
    for s in output_list:
        all_interactions = all_interactions.union(s)

    with open("predicted_ddi_ppi" + '.tsv', 'w') as f:
        for ddi in all_interactions:
            f.write(f'{ddi[0]}\t{ddi[1]}\n')
    print("Done")


if __name__ == '__main__':
    # write_ddi_ppi_connection()
    read_protein_interactions("../resultdata/source3_dippfam_ordered", inter_predicted, mart_table)
    # add_classification("predicted_ddi_ppi.tsv", "resultdata/result-all")
