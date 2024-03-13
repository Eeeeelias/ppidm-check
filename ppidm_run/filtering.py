import operator

import datetime
import pickle
import sys
import itertools

from ppidm_run.main import result_address, source_address
import random
from math import sqrt

# reproducibility yay
random.seed(42)


def read_interactions(file_path: str):
    interactions_3did = set()
    pfam_3did = set()
    file1 = open(file_path, 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            interactions_3did.add((item1, item2))
        else:
            interactions_3did.add((item2, item1))
        pfam_3did.add(item1)
        pfam_3did.add(item2)
    return interactions_3did, pfam_3did


def interactions(file_path: str, info: dict, info_tuple_multiple: dict, source: str) -> tuple[set, set, dict, dict]:
    interactions_intact = set()
    pfam_intact = set()
    file1 = open(file_path, 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            interaction = (line_sp[0], line_sp[1])
        else:
            interaction = (line_sp[1], line_sp[0])

        if '_' not in line:
            interactions_intact.add(interaction)
            if interaction[0] in info:
                if interaction[1] in info[interaction[0]]:
                    info[interaction[0]][interaction[1]][source] = float(line_sp[5])
                else:
                    info[interaction[0]][interaction[1]] = dict()
                    info[interaction[0]][interaction[1]][source] = float(line_sp[5])
            else:
                info[interaction[0]] = dict()
                info[interaction[0]][interaction[1]] = dict()
                info[interaction[0]][interaction[1]][source] = float(line_sp[5])

            pfam_intact.add(line_sp[0])
            pfam_intact.add(line_sp[1])
        else:
            if interaction[0] in info_tuple_multiple:
                if interaction[1] in info_tuple_multiple[interaction[0]]:
                    info_tuple_multiple[interaction[0]][interaction[1]][source] = float(line_sp[5])
                else:
                    info_tuple_multiple[interaction[0]][interaction[1]] = dict()
                    info_tuple_multiple[interaction[0]][interaction[1]][source] = float(line_sp[5])
            else:
                info_tuple_multiple[interaction[0]] = dict()
                info_tuple_multiple[interaction[0]][interaction[1]] = dict()
                info_tuple_multiple[interaction[0]][interaction[1]][source] = float(line_sp[5])
    return interactions_intact, pfam_intact, info, info_tuple_multiple


def coef_score(coefficients: list, interaction: list, sources: list):
    result = []
    for coef, src in zip(coefficients, sources):
        result.append(coef * interaction[src])
    return result


def extract_info(relevant_pfams: set, score_info: dict):
    result = {}
    for pfam_id in relevant_pfams:
        try:
            result[pfam_id] = score_info[pfam_id]
        except KeyError:
            result[pfam_id] = 0
    pickle.dump(result, open('../pickles/info_scores.pickle', 'wb'))
    print("Wrote scores to pickle")


def create_wrong_assocations(sources):
    start = datetime.datetime.now()
    print("Create NEGATIVE assocaitions from all inputs (around ? mins)")

    # removed kbdock from here

    file1 = open(result_address + '3did', 'r')
    gs = set()
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            interaction = (line_sp[0], line_sp[1])
        else:
            interaction = (line_sp[1], line_sp[0])
        gs.add(interaction)

    common_factors = set()
    dom_common_factors = dict()
    all_interactions_DDI = set()
    # for source in ['source1_intact']:
    for source in sources:
        file1 = open(result_address + source + 'pfam', 'r')
        # interaction_score = dict()
        for line in file1:
            line_sp = line.rstrip().split('\t')
            cf = hash(line_sp[1])
            dom1 = line_sp[0]
            dom2 = line_sp[2]
            if '_' in dom1 or '_' in dom2:
                continue

            all_interactions_DDI.add((dom1, dom2))

            if 'string' in source:
                continue

            common_factors.add(cf)

            if dom1 not in dom_common_factors:
                dom_common_factors[dom1] = set()
            dom_common_factors[dom1].add(cf)

            if dom2 not in dom_common_factors:
                dom_common_factors[dom2] = set()
            dom_common_factors[dom2].add(cf)

    wrongcf_foreach_dom = dict()
    print(len(common_factors))
    print(len(dom_common_factors))

    counter = 0
    used = set()
    result_file = open(result_address + 'negative_set', 'w')

    # for every domain that has interactions
    for dom in dom_common_factors:
        counter += 1
        # print(counter)
        # get the interactions of that domain
        number_cf_for_dom = dom_common_factors[dom]
        # remove them from all interactions in general
        all_possible_cf_for_dom = common_factors - number_cf_for_dom
        # randomly get the same amount of interactions from the rest (= node degree)
        wrong_cfs = random.sample(all_possible_cf_for_dom, len(number_cf_for_dom))
        # save it
        wrongcf_foreach_dom[dom] = wrong_cfs

    counter = 0
    for dom1 in dom_common_factors:
        for dom2 in dom_common_factors:
            counter += 1
            if counter % 100_000 == 0:
                print(counter)
            if dom1 < dom2:
                interaction = (dom1, dom2)
            elif dom1 > dom2:
                interaction = (dom2, dom1)
            else:
                continue

            if interaction in used:
                continue
            used.add(interaction)

            # making sure the interactions are not in the positive set already
            if interaction in gs or interaction in all_interactions_DDI:
                continue

            # set of all wrong common factors (hashes)
            wrong_cfs1: set = wrongcf_foreach_dom[dom1]

            wrong_cfs2: set = wrongcf_foreach_dom[dom2]

            intersection = set(wrong_cfs1) & set(wrong_cfs2)
            nominator = len(intersection)
            denom1 = len(set(wrong_cfs1))
            denom2 = len(set(wrong_cfs2))

            similarity = float(nominator) / (sqrt(denom1) * sqrt(denom2))
            if similarity != 0:
                result_file.write(f"{dom1}\t{dom2}\t{denom1}\t{denom2}\t{nominator}\t{similarity}\n")
    result_file.close()

    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")


def assign_interaction(sources):
    start = datetime.datetime.now()
    print("Filtering associations for Interactions(around ? mins)")
    interactions_3did, pfam_3did = read_interactions(result_address + '3did')
    # interactions_kbdock, pfam_kbdock = read_interactions(result_address + 'kbdock')
    # interactions_domine, pfam_domine = read_interactions(result_address + 'domine')

    info = dict()
    info_tuple_multiple = dict()

    interactions_sources = {}
    pfam_sources = {}
    source_names = []

    for source in sources:
        source_name = source[8:]
        print("Getting interactions for source:", source_name)
        interaction, pfam, info, info_tuple_multiple = interactions(
            result_address + 'pfam-pfam-interaction-' + source_name, info, info_tuple_multiple, source_name)
        source_names.append(source_name)
        interactions_sources[source_name] = interaction
        pfam_sources[source_name] = pfam

    for item1 in info:
        for item2 in info[item1]:
            for source in source_names:
                info[item1][item2].setdefault(source, 0)

    print("Getting gold_standard interactions: " + str(datetime.datetime.now() - start) + "\n")
    # find the best interactions that are well-supported by any of the sources
    gold_standard = set()
    for source in interactions_sources.keys():
        gold_standard = gold_standard.union(interactions_sources[source])
    background_data = list(gold_standard - (gold_standard.intersection(interactions_3did)))
    gold_standard = gold_standard.intersection(interactions_3did)

    train_set = random.sample(gold_standard, int(len(gold_standard) / 2))

    test_set = gold_standard - set(train_set)
    print("Length of gold standard:", len(gold_standard))
    print("Length of train set:", len(train_set))
    print("Length of test set:", len(test_set))
    print("Length of background data:", len(background_data))
    # backgroundData = random.sample(backgroundData, len(backgroundData) / 100)
    # print(len(backgroundData))

    best_area_under_curve = 0

    num_sources = len(sources)
    best_coefs = {f'coef{i}': 0 for i in range(1, num_sources + 1)}

    # TODO: give proper ranges that make sense for the coefficients and minimize the search space

    # coefficient_ranges = [range(1, 2), range(5, 6), range(1, 2), range(9, 10), range(12, 13), range(6, 7), range(100, 101),
    #                range(17, 18)]
    coefficient_ranges = [range(1, 2)] * num_sources

    # Generate all possible combinations of coefficients
    for coefficients in itertools.product(*coefficient_ranges):
        coefficients = list(coefficients)
        lenSize = len(background_data)
        coef_sum = sum(coefficients)
        all_data_positive_negative = {}
        for datum in gold_standard:
            if datum[0] in info:
                if datum[1] in info[datum[0]]:
                    dom_score = coef_score(coefficients, info[datum[0]][datum[1]], source_names)
                    all_data_positive_negative[datum] = sum(dom_score) / coef_sum
                else:
                    all_data_positive_negative[datum] = 0
            else:
                all_data_positive_negative[datum] = 0

        for datum in background_data:
            if datum[0] in info:
                if datum[1] in info[datum[0]]:
                    dom_score = coef_score(coefficients, info[datum[0]][datum[1]], source_names)
                    all_data_positive_negative[datum] = sum(dom_score) / coef_sum

        sorted_all_data_positive_negative = sorted(all_data_positive_negative.items(),
                                                   key=operator.itemgetter(1), reverse=True)
        yindex = 0.0
        area_under_curve = 0.0

        for datum in sorted_all_data_positive_negative:
            # print(datum[0])
            # print(gold_standard)
            # raw_input()
            if datum[0] in gold_standard:
                yindex += 1
            else:
                area_under_curve += (yindex / len(gold_standard)) * (1.0 / lenSize)

        print("AUC:", area_under_curve)
        print("Coefficients:", coefficients)
        print("Best Coefficients", list(best_coefs.values()))
        print("\n")

        if best_area_under_curve < area_under_curve:
            best_area_under_curve = area_under_curve
            best_coefs = {f'coef{i}': coef for i, coef in enumerate(coefficients, 1)}


    print("Overall best coefs:", list(best_coefs.values()))
    print("Best AUC:", best_area_under_curve)

    all_data_scores = dict()
    best_coefs = [x for x in best_coefs.values()]
    coef_summation = sum(best_coefs)
    for item1 in info:
        for item2 in info[item1]:
            dom_score = coef_score(best_coefs, info[item1][item2], source_names)
            if item1 in all_data_scores:
                all_data_scores[item1][item2] = sum(dom_score) / coef_summation
            else:
                all_data_scores[item1] = dict()
                all_data_scores[item1][item2] = sum(dom_score) / coef_summation

    # finding negative set with low-scoring
    print("finding low-scoring associations where score cap be anything")
    neg_model = 1

    gold_standard_negative_set = set()

    if neg_model == 1:
        negatives = set()
        negatives_score = dict()
        negative_file = open(result_address + 'negative_set', 'r')
        for line in negative_file:
            line_sp = line.rstrip().split('\t')
            # if float(line_sp[5]) > 0.04:
            #     continue
            negatives.add((line_sp[0], line_sp[1]))
            negatives_score[(line_sp[0], line_sp[1])] = float(line_sp[5])
            # print(line_sp[5])

        print("Negatives", len(negatives))
        gold_standard_negative_set = random.sample(negatives, len(gold_standard))

    elif neg_model == 2:
        max_pair = ''
        max_score = 0.0
        for datum1 in info:
            for datum2 in info[datum1]:
                if (datum1, datum2) in gold_standard:
                    continue
                flag_ok_for_negative = 0
                for source in source_names:
                    if info[datum1][datum2][source] > 0:
                        flag_ok_for_negative += 1

                if flag_ok_for_negative >= 6:
                    if len(gold_standard_negative_set) < len(gold_standard):
                        gold_standard_negative_set.add((datum1, datum2))
                        max_score = 0.0
                        for item in gold_standard_negative_set:
                            if all_data_scores[item[0]][item[1]] > max_score:
                                max_score = all_data_scores[item[0]][item[1]]
                                max_pair = item
                    else:
                        if all_data_scores[datum1][datum2] < max_score:
                            max_score = all_data_scores[datum1][datum2]
                            # print(max_pair)
                            gold_standard_negative_set.remove(max_pair)
                            gold_standard_negative_set.add((datum1, datum2))
                            max_score = 0.0
                            for item in gold_standard_negative_set:
                                if all_data_scores[item[0]][item[1]] > max_score:
                                    max_score = all_data_scores[item[0]][item[1]]
                                    max_pair = item
        print('negative set max score' + str(max_score))

    else:
        gold_standard_negative_set = random.sample(background_data, len(gold_standard))

    train_negative_set = random.sample(gold_standard_negative_set, int(len(gold_standard) / 2))
    test_negative_set = set(gold_standard_negative_set) - set(train_negative_set)

    print("Length of gold standard negative:", len(gold_standard_negative_set))
    print("Length of negative train:", len(train_negative_set))
    print("Length of negative test:", len(test_negative_set))
    # print(train_negative_set)

    best_fmeasure = 0
    best_threshold = 1000
    best_fmeasure_test = 0

    # plot distribution of scores for positives and negatives
    # calculating best Threshold and best F-measure
    count = 1
    for threshold in range(30, 1, -1):
        threshold = float(threshold) / 1000
        count_for_train = 0
        count_for_test = 0
        count_for_train_negative = 0
        count_for_test_negative = 0
        count_all_found = 0
        # Check the training set of interpro and see whether interproTrain association is found or not AND they are more
        # than THRESHOLD SCORE
        for datum in train_set:
            interaction = datum
            score = 0
            flag = False
            if datum[0] in all_data_scores:
                if datum[1] in all_data_scores[datum[0]]:

                    score = all_data_scores[datum[0]][datum[1]]
                    if score >= threshold:
                        # flag = True
                        count_for_train += 1

        for datum in test_set:
            interaction = datum
            score = 0
            # flag = False
            if datum[0] in all_data_scores:
                if datum[1] in all_data_scores[datum[0]]:

                    score = all_data_scores[datum[0]][datum[1]]
                    if score >= threshold:
                        # flag = True
                        count_for_test += 1

        for datum in train_negative_set:
            score = 0
            # flag = False
            if neg_model == 1:
                score = negatives_score[datum]
            else:
                if datum[0] in info:
                    if datum[1] in info[datum[0]]:
                        coef_summation = float(sum(best_coefs))
                        score = sum(coef_score(best_coefs, info[datum[0]][datum[1]], source_names)) / coef_summation
            if score >= threshold:
                # flag = True
                count_for_train_negative += 1


        for datum in test_negative_set:
            score = 0
            # flag = False
            if neg_model == 1:
                score = negatives_score[datum]
            else:
                if datum[0] in info:
                    if datum[1] in info[datum[0]]:
                        coef_summation = float(sum(best_coefs))
                        score = sum(coef_score(best_coefs, info[datum[0]][datum[1]], source_names)) / coef_summation
            if score >= threshold:
                # flag = True
                count_for_test_negative += 1

        # Check the found associations and see whether they are more than THRESHOLD SCORE
        for datum1 in info:
            for datum2 in info[datum1]:
                flag = False
                coef_summation = float(sum(best_coefs))
                score = sum(coef_score(best_coefs, info[datum1][datum2], source_names)) / coef_summation
                if score >= threshold:
                    flag = True
                    count_all_found += 1

        tp_train = count_for_train
        fn_train = len(train_set) - tp_train
        fp_train = count_for_train_negative
        f_score_train = (float(tp_train) * 2) / ((2 * tp_train) + fn_train + fp_train)

        tp_test = count_for_test
        fn_test = len(test_set) - tp_test
        fp_test = count_for_test_negative
        f_score_test = (float(tp_test) * 2) / ((2 * tp_test) + fn_test + fp_test)

        if best_fmeasure < f_score_train:
            best_fmeasure = f_score_train
            best_threshold = threshold
            best_fmeasure_test = f_score_test

        bold = '\033[1m'
        end = '\033[0m'
        print(f"############### Iteration {count} #################")
        print('Training Set'.ljust(50), "| Testing Set")
        print(f"T\t\tPredicted".ljust(45), "| T\t\tPredicted")
        print("T\t\tPos\t   Neg".ljust(45), "| T\t\tPos\t   Neg")
        print(f"T Pos\t{tp_train:4} | {fn_train}".ljust(48), f"| T Pos\t{tp_test} | {fn_test}")
        print(f"T Neg\t{fp_train:4} | / ".ljust(48), f"| T Neg\t{fp_test:4} | / ")
        print(f"T F_score: {round(f_score_train, 5)} | Best F_measure: {bold + str(round(best_fmeasure, 5)) + end}".ljust(58),
              f"| T F_score: {round(f_score_test, 5)} | Best F_measure: {bold + str(round(best_fmeasure_test, 5)) + end}")
        print("-" * 80)
        print(
            f"Current threshold: {threshold} | Best threshold: {bold + str(best_threshold) + end} | Count all found: {count_all_found}\n")
        count += 1

    result_calculated = open(result_address + 'pfam-pfam-interaction-calculated', 'w')
    result_merged = open(result_address + 'pfam-pfam-interaction-merged', 'w')
    for datum1 in info:
        for datum2 in info[datum1]:
            score = all_data_scores[datum1][datum2]
            source_infos = '\t'.join([str(info[datum1][datum2][source]) for source in source_names])
            result_string = datum1 + '\t' + datum2 + '\t' + source_infos + '\t' + str(score) + '\n'
            result_merged.write(result_string)
            if score >= best_threshold:
                flag = True
                result_calculated.write(result_string)

    result_gold_standard = open(result_address + 'pfam-pfam-interaction-goldstandard', 'w')
    for datum in gold_standard:
        if datum in train_set:
            flag = 'train-'
        else:
            flag = 'test-'

        try:
            score = all_data_scores[datum[0]][datum[1]]
        except KeyError:
            score = 0

        if score >= best_threshold:
            flag = flag + 'yes'
        else:
            flag = flag + 'no'

        sources = '\t'.join([str(info[datum[0]][datum[1]][source]) for source in source_names])
        result_gold = datum[0] + '\t' + datum[1] + '\t' + sources + '\t' + str(score) + '\t' + flag + '\n'
        result_gold_standard.write(result_gold)

    result_negative = open(result_address + 'pfam-pfam-interaction-negative', 'w')
    for datum in gold_standard_negative_set:
        if datum in train_negative_set:
            flag = 'train-'
        else:
            flag = 'test-'

        if neg_model == 1:
            score = negatives_score[datum]
        else:
            coef_summation = sum(best_coefs)
            score = sum(coef_score(best_coefs, info[datum[0]][datum[1]], source_names)) / coef_summation

        if score >= best_threshold:
            flag = flag + 'yes'
        else:
            flag = flag + 'no'

        result_negative.write(
            datum[0] + '\t' + datum[1] + '\t' + str(score) + '\t' + flag + '\n')

    for item1 in info_tuple_multiple:
        for item2 in info_tuple_multiple[item1]:
            for source in source_names:
                info_tuple_multiple[item1][item2].setdefault(source, 0)

    result_calculated_tuple = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'w')
    result_merged_tuple = open(result_address + 'pfam-pfam-interaction-merged_tuple', 'w')
    for datum1 in info_tuple_multiple:
        for datum2 in info_tuple_multiple[datum1]:
            # flag = False
            coef_summation = sum(best_coefs)
            score = sum(coef_score(best_coefs, info_tuple_multiple[datum1][datum2], source_names)) / coef_summation

            source_infos = '\t'.join([str(info_tuple_multiple[datum1][datum2][source]) for source in source_names])
            result_string = datum1 + '\t' + datum2 + '\t' + source_infos + '\t' + str(score) + '\n'

            result_merged_tuple.write(result_string)

            if score >= best_threshold:
                flag = True
                result_calculated_tuple.write(result_string)

    result_calculated.close()
    result_merged.close()
    result_calculated_tuple.close()
    result_merged_tuple.close()
    result_gold_standard.close()
    result_negative.close()
    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")
