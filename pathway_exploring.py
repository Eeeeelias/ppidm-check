import pickle
import pandas as pd


def read_pathways(file: str):
    pathway = []
    with open(file, 'r') as f:
        for line in f.readlines():
            if line[2:5] != 'HSA':
                continue
            else:
                split = line.strip().split("\t")
                pathway.append((split[0], split[1]))
    # print(len(pathway))
    return pathway


def relevant_subpathways(top_pathway: str, hierarchy: list[tuple[str, str]]):
    return [bot for top, bot in hierarchy if top == top_pathway]


def identified_pathways(enriched):
    top_enriched = enriched[enriched['Pathway ID'].isin(head_pathways)]

    print("High level enriched pathways and their sub-pathways:")
    for val in top_enriched['Pathway ID'].values:
        subpath = relevant_subpathways(val, pathway_hierarchy)
        print(f"{pathway_names[val]}\n", enriched[enriched['Pathway ID'].isin(subpath)]['Pathway name'].to_string(),
              sep="")


pathway_hierarchy = read_pathways('data_nease/ReactomePathwaysRelation.txt')
pathway_names = {key: val for key, val in read_pathways('data_nease/ReactomePathways.txt')}

head_pathways = set([a[0] for a in pathway_hierarchy]) - set([a[1] for a in pathway_hierarchy])
print("Top pathways:", len(head_pathways))
lowest_pathways = set([a[1] for a in pathway_hierarchy]) - set([a[0] for a in pathway_hierarchy])
print("Bottom pathways:", len(lowest_pathways))

enr_ext: pd.DataFrame = pickle.load(open('data_nease/nease_enr_ext.pkl', 'rb'))
enr_org: pd.DataFrame = pickle.load(open('data_nease/nease_enr_org.pkl', 'rb'))
cutoff = 0.05

enr_ext = enr_ext.loc[enr_ext['adj p_value'] <= cutoff]
enr_org = enr_org.loc[enr_org['adj p_value'] <= cutoff]

print("Enriched pathways in extended:", len(enr_ext['Pathway ID']))
print("Enriched pathways in original:", len(enr_org['Pathway ID']))

original_only = set(enr_org['Pathway ID']) - set(enr_ext['Pathway ID'])
print("Enriched pathways exclusively in original:", len(original_only), "\n")

identified_pathways(enr_ext)

bot_enriched = enr_ext[enr_ext['Pathway ID'].isin(lowest_pathways)]
print(f"\nHighest enriched lowest level pathways:\n#  {'Pathway name':63} {'Pathway ID':13} {'Nease score'}")
for idx, val in bot_enriched.head(10).iterrows():
    print(f"{idx:2} {val['Pathway ID']:13} {val['Pathway name']:63} {round(val['Nease score'], 3)}")
