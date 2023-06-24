import pickle
import seaborn as sns
from matplotlib import pyplot as plt

scores = pickle.load(open('pickles/info_scores.pickle', 'rb'))
gold_interactions = pickle.load(open('pickles/train_set.pickle', 'rb'))
print("Gold interactions:", len(gold_interactions))
random_interactions = pickle.load(open('pickles/random_train.pickle', 'rb'))
print("Random interactions:", len(random_interactions))
negative_interactions = pickle.load(open('pickles/negative_train.pickle', 'rb'))
print("Negative train interactions:", len(negative_interactions))


def get_scores(source, interactions):
    found_scores = []
    for interaction in source:
        if interaction[0] in interactions:
            if interactions[interaction[0]] == 0:
                found_scores.append(0)
                continue
            if interaction[1] in interactions[interaction[0]]:
                found_scores.append(interactions[interaction[0]][interaction[1]])
            else:
                found_scores.append(0)
        else:
            found_scores.append(0)
    # print(len(found_scores))
    return found_scores


gold_scores = get_scores(gold_interactions, scores)
random_scores = get_scores(random_interactions, scores)
negative_scores = list(negative_interactions.values())

# plotting
plt.style.use('ggplot')
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

sns.kdeplot(negative_scores, color='orange', ax=ax1)
ax1.axvline(x=0.01586, linestyle='dotted', color='#0097ff', label='threshold (0.01586)')
ax1.axvline(x=0.03, linestyle='dotted', color='#fc1b80', label='threshold (0.03)')
ax1.set_ylabel('Negative train - Density')


sns.kdeplot(gold_scores, color='blue', ax=ax2)
ax2.axvline(x=0.01586, linestyle='dotted', color='#0097ff', label='threshold')
ax2.axvline(x=0.03, linestyle='dotted', color='#fc1b80')
ax2.set_ylabel('Positive train - Density')

plt.xlim(0, 1)

# Adding labels to the plot
plt.xlabel('Score')
plt.suptitle('Score density in negatives vs. positives')

# Displaying the plot
plt.tight_layout()
ax1.legend()
plt.savefig('pictures/score_density.png')
plt.show()
