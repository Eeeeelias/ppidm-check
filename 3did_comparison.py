
def read_interactions(file: str):
    interactions = []
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.split("\t")
            interactions.append((line[0], line[1]))
    return interactions


did_2022 = read_interactions('resultdata/3did_2022')
did_2017 = read_interactions('resultdata/3did')
interactions_gold = read_interactions('resultdata/interactions_gold')

print("Length of interactions_gold:", len(interactions_gold))
print("Length of did_2022:", len(did_2022))
print("Length of did_2017:", len(did_2017))
print("Unique 2022:", len(set(did_2017) - set(did_2022)))
print("Interaction overlap gold and did_2017:", len(set(interactions_gold) & set(did_2017)))
print("Interaction overlap gold and did_2022:", len(set(interactions_gold) & set(did_2022)))
