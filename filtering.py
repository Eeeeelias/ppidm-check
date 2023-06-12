import operator

import datetime
import pickle

from main import result_address, source_address
import random
from math import sqrt


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


def create_wrong_assocations():
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
    for source in ['source1_intact', 'source2_mint', 'source3_dip', 'source4_biogrid', 'source5_string-exp',
                   'source5_string-rest', 'source6_sifts', 'source7_hprd']:
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

            if dom1 in dom_common_factors:
                dom_common_factors[dom1].add(cf)
            else:
                dom_common_factors[dom1] = set()
            dom_common_factors[dom1].add(cf)

            if dom2 in dom_common_factors:
                dom_common_factors[dom2].add(cf)
            else:
                dom_common_factors[dom2] = set()
            dom_common_factors[dom2].add(cf)

    wrongcf_foreach_dom = dict()
    print(len(common_factors))
    print(len(dom_common_factors))

    counter = 0
    used = set()
    result_file = open(result_address + 'negative_set1', 'w')

    for dom in dom_common_factors:
        counter += 1
        # print(counter)
        number_cf_for_dom = dom_common_factors[dom]
        all_possible_cf_for_dom = common_factors - number_cf_for_dom
        wrong_cfs = random.sample(all_possible_cf_for_dom, len(number_cf_for_dom))
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

            if interaction in gs or interaction in all_interactions_DDI:
                continue

            # number_cf_for_dom1 = dom_common_factors[dom1]
            # all_possible_cf_for_dom1 = common_factors - number_cf_for_dom1
            wrong_cfs1 = wrongcf_foreach_dom[dom1]

            # number_cf_for_dom2 = dom_common_factors[dom2]
            # all_possible_cf_for_dom2 = common_factors - number_cf_for_dom2
            wrong_cfs2 = wrongcf_foreach_dom[dom2]

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


def assign_interaction():
    start = datetime.datetime.now()
    print("Filtering associations for Interactions(around ? mins)")
    interactions_3did, pfam_3did = read_interactions(result_address + '3did')
    # interactions_kbdock, pfam_kbdock = read_interactions(result_address + 'kbdock')
    # interactions_domine, pfam_domine = read_interactions(result_address + 'domine')

    info = dict()
    info_tuple_multiple = dict()

    interactions_intact, pfam_intact, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-intact', info, info_tuple_multiple, 'intact')

    interactions_dip, pfam_dip, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-dip', info, info_tuple_multiple, 'dip')

    interactions_mint, pfam_mint, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-mint', info, info_tuple_multiple, 'mint')

    interactions_biogrid, pfam_biogrid, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-biogrid', info, info_tuple_multiple, 'biogrid')

    interactions_stringg_exp, pfam_stringg_exp, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-string-exp', info, info_tuple_multiple, 'string_exp')

    interactions_stringg_rest, pfam_stringg_rest, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-string-rest', info, info_tuple_multiple, 'string_rest')

    interactions_sifts_accession, pfam_sifts_accession, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-sifts', info, info_tuple_multiple, 'sifts_acc')

    interactions_hprd, pfam_hprd, info, info_tuple_multiple = interactions(
        result_address + 'pfam-pfam-interaction-hprd', info, info_tuple_multiple, 'hprd')

    for item1 in info:
        for item2 in info[item1]:
            if 'intact' not in info[item1][item2]:
                info[item1][item2]['intact'] = 0
            if 'dip' not in info[item1][item2]:
                info[item1][item2]['dip'] = 0
            if 'mint' not in info[item1][item2]:
                info[item1][item2]['mint'] = 0
            if 'biogrid' not in info[item1][item2]:
                info[item1][item2]['biogrid'] = 0
            if 'string_exp' not in info[item1][item2]:
                info[item1][item2]['string_exp'] = 0
            if 'string_rest' not in info[item1][item2]:
                info[item1][item2]['string_rest'] = 0
            if 'sifts_acc' not in info[item1][item2]:
                info[item1][item2]['sifts_acc'] = 0
            if 'hprd' not in info[item1][item2]:
                info[item1][item2]['hprd'] = 0

    # removed interactions_kbdock here due to me not having these interactions
    gold_standard = (interactions_intact |
                     interactions_dip |
                     interactions_mint |
                     interactions_biogrid |
                     interactions_stringg_exp |
                     interactions_stringg_rest |
                     interactions_sifts_accession |
                     interactions_hprd) & \
                    (interactions_3did)

    print("Getting gold_standard interactions: " + str(datetime.datetime.now() - start) + "\n")
    pickle.dump(gold_standard, open('gold_standard.pickle', 'wb'))

    train_set = random.sample(gold_standard, int(len(gold_standard) / 2))
    test_set = gold_standard - set(train_set)
    print("Length of gold standard:", len(gold_standard))
    print("Length of train set:", len(train_set))
    print("Length of test set:", len(test_set))

    backgroundData = []
    for item1 in info:
        for item2 in info[item1]:
            if (item1, item2) not in gold_standard:
                backgroundData.append((item1, item2))

    print(len(backgroundData))
    # backgroundData = random.sample(backgroundData, len(backgroundData) / 100)
    # print(len(backgroundData))

    best_area_under_curve = 0
    best_coef1 = 0
    best_coef2 = 0
    best_coef3 = 0
    best_coef4 = 0
    best_coef5 = 0
    best_coef6 = 0
    best_coef7 = 0
    best_coef8 = 0

    for coef1 in range(5, 6):
        for coef2 in range(1, 2):
            for coef3 in range(1, 2):
                for coef4 in range(9, 10):
                    for coef5 in range(12, 13):
                        for coef6 in range(6, 7):
                            for coef7 in range(100, 101):
                                for coef8 in range(17, 18):
                                    lenSize = len(backgroundData)
                                    coef_summation = coef1 + coef2 + coef3 + coef4 + coef5 + coef6 + coef7 + coef8
                                    all_data_positive_negative = {}
                                    for datum in gold_standard:
                                        if datum[0] in info:
                                            if datum[1] in info[datum[0]]:
                                                first_part = coef1 * info[datum[0]][datum[1]]['intact']
                                                second_part = coef2 * info[datum[0]][datum[1]]['dip']
                                                third_part = coef3 * info[datum[0]][datum[1]]['mint']
                                                fourth_part = coef4 * info[datum[0]][datum[1]]['biogrid']
                                                fifth_part = coef5 * info[datum[0]][datum[1]]['string_exp']
                                                sixth_part = coef6 * info[datum[0]][datum[1]]['string_rest']
                                                seventh_part = coef7 * info[datum[0]][datum[1]]['sifts_acc']
                                                eighth_part = coef8 * info[datum[0]][datum[1]]['hprd']
                                                all_data_positive_negative[datum] = (
                                                                                            first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
                                            else:
                                                all_data_positive_negative[datum] = 0
                                        else:
                                            all_data_positive_negative[datum] = 0

                                    for datum in backgroundData:
                                        if datum[0] in info:
                                            if datum[1] in info[datum[0]]:
                                                first_part = coef1 * info[datum[0]][datum[1]]['intact']
                                                second_part = coef2 * info[datum[0]][datum[1]]['dip']
                                                third_part = coef3 * info[datum[0]][datum[1]]['mint']
                                                fourth_part = coef4 * info[datum[0]][datum[1]]['biogrid']
                                                fifth_part = coef5 * info[datum[0]][datum[1]]['string_exp']
                                                sixth_part = coef6 * info[datum[0]][datum[1]]['string_rest']
                                                seventh_part = coef7 * info[datum[0]][datum[1]]['sifts_acc']
                                                eighth_part = coef8 * info[datum[0]][datum[1]]['hprd']
                                                all_data_positive_negative[datum] = (
                                                                                            first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation

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

                                    print(area_under_curve)
                                    print(coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8)
                                    print(
                                        best_coef1, best_coef2, best_coef3, best_coef4, best_coef5, best_coef6,
                                        best_coef7, best_coef8)
                                    print("\n")

                                    if best_area_under_curve < area_under_curve:
                                        best_area_under_curve = area_under_curve
                                        best_coef1 = coef1
                                        best_coef2 = coef2
                                        best_coef3 = coef3
                                        best_coef4 = coef4
                                        best_coef5 = coef5
                                        best_coef6 = coef6
                                        best_coef7 = coef7
                                        best_coef8 = coef8

    print(best_coef1, best_coef2, best_coef3, best_coef4, best_coef5, best_coef6, best_coef7, best_coef8)
    print(best_area_under_curve)

    all_data_scores = dict()
    coef_summation = best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8
    for item1 in info:
        for item2 in info[item1]:
            first_part = best_coef1 * info[item1][item2]['intact']
            second_part = best_coef2 * info[item1][item2]['dip']
            third_part = best_coef3 * info[item1][item2]['mint']
            fourth_part = best_coef4 * info[item1][item2]['biogrid']
            fifth_part = best_coef5 * info[item1][item2]['string_exp']
            sixth_part = best_coef6 * info[item1][item2]['string_rest']
            seventh_part = best_coef7 * info[item1][item2]['sifts_acc']
            eighth_part = best_coef8 * info[item1][item2]['hprd']
            if item1 in all_data_scores:
                all_data_scores[item1][item2] = (
                                                        first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
            else:
                all_data_scores[item1] = dict()
                all_data_scores[item1][item2] = (
                                                        first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation

    # finding negative set with low-scoring
    print("finding low-scoring associations where score cap be anything")
    neg_model = 1

    gold_standard_negative_set = set()

    if neg_model == 1:
        negatives = set()
        negatives_score = dict()
        negative_file = open(result_address + 'negative_set1', 'r')
        for line in negative_file:
            line_sp = line.rstrip().split('\t')
            # if float(line_sp[5]) > 0.04:
            #     continue
            negatives.add((line_sp[0], line_sp[1]))
            negatives_score[(line_sp[0], line_sp[1])] = float(line_sp[5])
            # print(line_sp[5])

        print(len(negatives))
        gold_standard_negative_set = random.sample(negatives, len(gold_standard))

    elif neg_model == 2:
        max_pair = ''
        max_score = 0.0
        for datum1 in info:
            for datum2 in info[datum1]:
                if (datum1, datum2) in gold_standard:
                    continue
                flag_ok_for_negative = 0
                if info[datum1][datum2]['intact'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['dip'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['mint'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['biogrid'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['string_exp'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['string_rest'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['sifts_acc'] > 0:
                    flag_ok_for_negative += 1
                if info[datum1][datum2]['hprd'] > 0:
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
        gold_standard_negative_set = random.sample(backgroundData, len(gold_standard))

    train_negative_set = random.sample(gold_standard_negative_set, int(len(gold_standard) / 2))
    test_negative_set = set(gold_standard_negative_set) - set(train_negative_set)

    print("Length of gold standard negative:", len(gold_standard_negative_set))
    print("Length of negative train:", len(train_negative_set))
    print("Length of negative test:", len(test_negative_set))

    # print(train_negative_set)

    best_fmeasure = 0
    best_threshold = 1000
    best_fmeasure_test = 0

    # calculating best Threshold and best F-measure
    for threshold in range(30, 1, -1):
        threshold = float(threshold) / 1000
        count_for_train = 0
        count_for_test = 0
        count_for_train_negative = 0
        count_for_test_negative = 0
        count_all_found = 0
        # Check the training set of interpro and see whether interproTrain association is found or not AND they are more than THRESHOLD SCORE
        for datum in train_set:
            interaction = datum
            score = 0
            flag = False
            if datum[0] in all_data_scores:
                if datum[1] in all_data_scores[datum[0]]:
                    # coef_summation = best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8
                    # first_part = best_coef1 * info[datum[0]][datum[1]]['intact']
                    # second_part = best_coef2 * info[datum[0]][datum[1]]['dip']
                    # third_part = best_coef3 * info[datum[0]][datum[1]]['mint']
                    # fourth_part = best_coef4 * info[datum[0]][datum[1]]['biogrid']
                    # fifth_part = best_coef5 * info[datum[0]][datum[1]]['string_exp']
                    # sixth_part = best_coef6 * info[datum[0]][datum[1]]['string_rest']
                    # seventh_part = best_coef7 * info[datum[0]][datum[1]]['sifts_acc']
                    # eighth_part = best_coef8 * info[datum[0]][datum[1]]['hprd']
                    score = all_data_scores[datum[0]][datum[1]]
                    # score = (
                    #             first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
                    if score >= threshold:
                        # flag = True
                        count_for_train += 1

        for datum in test_set:
            interaction = datum
            score = 0
            # flag = False
            if datum[0] in all_data_scores:
                if datum[1] in all_data_scores[datum[0]]:
                    # coef_summation = float(best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8)
                    # first_part = best_coef1 * info[datum[0]][datum[1]]['intact']
                    # second_part = best_coef2 * info[datum[0]][datum[1]]['dip']
                    # third_part = best_coef3 * info[datum[0]][datum[1]]['mint']
                    # fourth_part = best_coef4 * info[datum[0]][datum[1]]['biogrid']
                    # fifth_part = best_coef5 * info[datum[0]][datum[1]]['string_exp']
                    # sixth_part = best_coef6 * info[datum[0]][datum[1]]['string_rest']
                    # seventh_part = best_coef7 * info[datum[0]][datum[1]]['sifts_acc']
                    # eighth_part = best_coef8 * info[datum[0]][datum[1]]['hprd']

                    score = all_data_scores[datum[0]][datum[1]]
                    # score = (
                    #             first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
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
                        coef_summation = float(
                            best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8)
                        first_part = best_coef1 * info[datum[0]][datum[1]]['intact']
                        second_part = best_coef2 * info[datum[0]][datum[1]]['dip']
                        third_part = best_coef3 * info[datum[0]][datum[1]]['mint']
                        fourth_part = best_coef4 * info[datum[0]][datum[1]]['biogrid']
                        fifth_part = best_coef5 * info[datum[0]][datum[1]]['string_exp']
                        sixth_part = best_coef6 * info[datum[0]][datum[1]]['string_rest']
                        seventh_part = best_coef7 * info[datum[0]][datum[1]]['sifts_acc']
                        eighth_part = best_coef8 * info[datum[0]][datum[1]]['hprd']

                        score = (
                                        first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
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
                        coef_summation = float(
                            best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8)
                        first_part = best_coef1 * info[datum[0]][datum[1]]['intact']
                        second_part = best_coef2 * info[datum[0]][datum[1]]['dip']
                        third_part = best_coef3 * info[datum[0]][datum[1]]['mint']
                        fourth_part = best_coef4 * info[datum[0]][datum[1]]['biogrid']
                        fifth_part = best_coef5 * info[datum[0]][datum[1]]['string_exp']
                        sixth_part = best_coef6 * info[datum[0]][datum[1]]['string_rest']
                        seventh_part = best_coef7 * info[datum[0]][datum[1]]['sifts_acc']
                        eighth_part = best_coef8 * info[datum[0]][datum[1]]['hprd']

                        score = (
                                        first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
            if score >= threshold:
                # flag = True
                count_for_test_negative += 1

        # Check the found associations and see whether they are more than THRESHOLD SCORE
        for datum1 in info:
            for datum2 in info[datum1]:
                flag = False
                coef_summation = float(
                    best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8)
                first_part = best_coef1 * info[datum1][datum2]['intact']
                second_part = best_coef2 * info[datum1][datum2]['dip']
                third_part = best_coef3 * info[datum1][datum2]['mint']
                fourth_part = best_coef4 * info[datum1][datum2]['biogrid']
                fifth_part = best_coef5 * info[datum1][datum2]['string_exp']
                sixth_part = best_coef6 * info[datum1][datum2]['string_rest']
                seventh_part = best_coef7 * info[datum1][datum2]['sifts_acc']
                eighth_part = best_coef8 * info[datum1][datum2]['hprd']

                score = (
                                first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
                if score >= threshold:
                    flag = True
                    count_all_found += 1

        # print(count_all_found)
        # print(count_for_train)
        # print(count_for_test)
        # print(count_for_train_negative)
        # print(count_for_test_negative)
        # print(' ')

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

        print(tp_train, fn_train, fp_train, f_score_train, best_fmeasure)
        print(threshold, best_threshold)
        print(tp_test, fn_test, fp_test, f_score_test, best_fmeasure_test)
        print(count_all_found)
        print(' ')

    result_calculated = open(result_address + 'pfam-pfam-interaction-calculated', 'w')
    result_merged = open(result_address + 'pfam-pfam-interaction-merged', 'w')
    for datum1 in info:
        for datum2 in info[datum1]:
            # flag = False
            # coef_summation = best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8
            # first_part = best_coef1 * info[datum1][datum2]['intact']
            # second_part = best_coef2 * info[datum1][datum2]['dip']
            # third_part = best_coef3 * info[datum1][datum2]['mint']
            # fourth_part = best_coef4 * info[datum1][datum2]['biogrid']
            # fifth_part = best_coef5 * info[datum1][datum2]['string_exp']
            # sixth_part = best_coef6 * info[datum1][datum2]['string_rest']
            # seventh_part = best_coef7 * info[datum1][datum2]['sifts_acc']
            # eighth_part = best_coef8 * info[datum1][datum2]['hprd']
            score = all_data_scores[datum1][datum2]
            # score = (
            #             first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
            result_merged.write(datum1 + '\t' + datum2 + '\t' + str(info[datum1][datum2]['intact']) + '\t' + str(
                info[datum1][datum2]['dip']) + '\t' + str(info[datum1][datum2]['mint']) + '\t' + str(
                info[datum1][datum2]['biogrid']) + '\t' + str(info[datum1][datum2]['string_exp']) + '\t' + str(
                info[datum1][datum2]['string_rest']) + '\t' + str(info[datum1][datum2]['sifts_acc']) + '\t' + str(
                info[datum1][datum2]['hprd']) + '\t' + str(
                score) + '\n')

            if score >= best_threshold:
                flag = True
                result_calculated.write(
                    datum1 + '\t' + datum2 + '\t' + str(info[datum1][datum2]['intact']) + '\t' + str(
                        info[datum1][datum2]['dip']) + '\t' + str(info[datum1][datum2]['mint']) + '\t' + str(
                        info[datum1][datum2]['biogrid']) + '\t' + str(info[datum1][datum2]['string_exp']) + '\t' + str(
                        info[datum1][datum2]['string_rest']) + '\t' + str(
                        info[datum1][datum2]['sifts_acc']) + '\t' + str(
                        info[datum1][datum2]['hprd']) + '\t' + str(
                        score) + '\n')

    result_gold_standard = open(result_address + 'pfam-pfam-interaction-goldstandard', 'w')
    for datum in gold_standard:
        if datum in train_set:
            flag = 'train-'
        else:
            flag = 'test-'

        # coef_summation = best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8
        # first_part = best_coef1 * info[datum[0]][datum[1]]['intact']
        # second_part = best_coef2 * info[datum[0]][datum[1]]['dip']
        # third_part = best_coef3 * info[datum[0]][datum[1]]['mint']
        # fourth_part = best_coef4 * info[datum[0]][datum[1]]['biogrid']
        # fifth_part = best_coef5 * info[datum[0]][datum[1]]['string_exp']
        # sixth_part = best_coef6 * info[datum[0]][datum[1]]['string_rest']
        # seventh_part = best_coef7 * info[datum[0]][datum[1]]['sifts_acc']
        # eighth_part = best_coef8 * info[datum[0]][datum[1]]['hprd']
        score = all_data_scores[datum[0]][datum[1]]
        # score = (
        #             first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation

        if score >= best_threshold:
            flag = flag + 'yes'
        else:
            flag = flag + 'no'

        result_gold_standard.write(
            datum[0] + '\t' + datum[1] + '\t' + str(info[datum[0]][datum[1]]['intact']) + '\t' + str(
                info[datum[0]][datum[1]]['dip']) + '\t' + str(info[datum[0]][datum[1]]['mint']) + '\t' + str(
                info[datum[0]][datum[1]]['biogrid']) + '\t' + str(info[datum[0]][datum[1]]['string_exp']) + '\t' + str(
                info[datum[0]][datum[1]]['string_rest']) + '\t' + str(
                info[datum[0]][datum[1]]['sifts_acc']) + '\t' + str(info[datum[0]][datum[1]]['hprd']) + '\t' + str(
                score) + '\t' + flag + '\n')

    result_negative = open(result_address + 'pfam-pfam-interaction-negative', 'w')
    for datum in gold_standard_negative_set:
        if datum in train_negative_set:
            flag = 'train-'
        else:
            flag = 'test-'

        if neg_model == 1:
            score = negatives_score[datum]
        else:
            coef_summation = best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8
            first_part = best_coef1 * info[datum[0]][datum[1]]['intact']
            second_part = best_coef2 * info[datum[0]][datum[1]]['dip']
            third_part = best_coef3 * info[datum[0]][datum[1]]['mint']
            fourth_part = best_coef4 * info[datum[0]][datum[1]]['biogrid']
            fifth_part = best_coef5 * info[datum[0]][datum[1]]['string_exp']
            sixth_part = best_coef6 * info[datum[0]][datum[1]]['string_rest']
            seventh_part = best_coef7 * info[datum[0]][datum[1]]['sifts_acc']
            eighth_part = best_coef8 * info[datum[0]][datum[1]]['hprd']

            score = (
                            first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation

        if score >= best_threshold:
            flag = flag + 'yes'
        else:
            flag = flag + 'no'

        result_negative.write(
            datum[0] + '\t' + datum[1] + '\t' + str(score) + '\t' + flag + '\n')

    for item1 in info_tuple_multiple:
        for item2 in info_tuple_multiple[item1]:
            if 'intact' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['intact'] = 0
            if 'dip' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['dip'] = 0
            if 'mint' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['mint'] = 0
            if 'biogrid' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['biogrid'] = 0
            if 'string_exp' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['string_exp'] = 0
            if 'string_rest' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['string_rest'] = 0
            if 'sifts_acc' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['sifts_acc'] = 0
            if 'hprd' not in info_tuple_multiple[item1][item2]:
                info_tuple_multiple[item1][item2]['hprd'] = 0

    result_calculated_tuple = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'w')
    result_merged_tuple = open(result_address + 'pfam-pfam-interaction-merged_tuple', 'w')
    for datum1 in info_tuple_multiple:
        for datum2 in info_tuple_multiple[datum1]:
            # flag = False
            coef_summation = best_coef1 + best_coef2 + best_coef3 + best_coef4 + best_coef5 + best_coef6 + best_coef7 + best_coef8
            first_part = best_coef1 * info_tuple_multiple[datum1][datum2]['intact']
            second_part = best_coef2 * info_tuple_multiple[datum1][datum2]['dip']
            third_part = best_coef3 * info_tuple_multiple[datum1][datum2]['mint']
            fourth_part = best_coef4 * info_tuple_multiple[datum1][datum2]['biogrid']
            fifth_part = best_coef5 * info_tuple_multiple[datum1][datum2]['string_exp']
            sixth_part = best_coef6 * info_tuple_multiple[datum1][datum2]['string_rest']
            seventh_part = best_coef7 * info_tuple_multiple[datum1][datum2]['sifts_acc']
            eighth_part = best_coef8 * info_tuple_multiple[datum1][datum2]['hprd']
            score = (
                            first_part + second_part + third_part + fourth_part + fifth_part + sixth_part + seventh_part + eighth_part) / coef_summation
            result_merged_tuple.write(
                datum1 + '\t' + datum2 + '\t' + str(info_tuple_multiple[datum1][datum2]['intact']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['dip']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['mint']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['biogrid']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['string_exp']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['string_rest']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['sifts_acc']) + '\t' + str(
                    info_tuple_multiple[datum1][datum2]['hprd']) + '\t' + str(
                    score) + '\n')

            if score >= best_threshold:
                flag = True
                result_calculated_tuple.write(
                    datum1 + '\t' + datum2 + '\t' + str(info_tuple_multiple[datum1][datum2]['intact']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['dip']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['mint']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['biogrid']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['string_exp']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['string_rest']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['sifts_acc']) + '\t' + str(
                        info_tuple_multiple[datum1][datum2]['hprd']) + '\t' + str(
                        score) + '\n')

    result_calculated.close()
    result_merged.close()
    result_calculated_tuple.close()
    result_merged_tuple.close()
    result_gold_standard.close()
    result_negative.close()
    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")
