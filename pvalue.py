import datetime
import math
from math import sqrt

from main import source_address, result_address


def pvalue_calculation(source, seqDom, pdbchainDom):
    start = datetime.datetime.now()
    print("P-value calculation for %s (around 2 secs)" % source)

    allPPI = set()
    file1 = open(source_address + source, 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            seq_seq = item1 + '_' + item2
        else:
            seq_seq = item2 + '_' + item1

        allPPI.add(seq_seq)

    N = len(allPPI)

    print(N)

    domain_seq_seq = dict()
    for seq_seq in allPPI:
        seq_split = seq_seq.split('_')
        if len(seq_split) > 2:
            seq1 = seq_split[0] + '_' + seq_split[1]
            seq2 = seq_split[2] + '_' + seq_split[3]
        else:
            seq1 = seq_split[0]
            seq2 = seq_split[1]
        if seq1 in seqDom:
            domains = seqDom[seq1]
            for domain in domains:
                if domain in domain_seq_seq:
                    domain_seq_seq[domain].add(seq_seq)
                else:
                    domain_seq_seq[domain] = set()
                    domain_seq_seq[domain].add(seq_seq)

        if seq2 in seqDom:
            domains = seqDom[seq2]
            for domain in domains:
                if domain in domain_seq_seq:
                    domain_seq_seq[domain].add(seq_seq)
                else:
                    domain_seq_seq[domain] = set()
                    domain_seq_seq[domain].add(seq_seq)

        if source == 'source6_sifts':
            if seq1 in pdbchainDom:
                domains = pdbchainDom[seq1]
                for domain in domains:
                    if domain in domain_seq_seq:
                        domain_seq_seq[domain].add(seq_seq)
                    else:
                        domain_seq_seq[domain] = set()
                        domain_seq_seq[domain].add(seq_seq)

            if seq2 in pdbchainDom:
                domains = pdbchainDom[seq2]
                for domain in domains:
                    if domain in domain_seq_seq:
                        domain_seq_seq[domain].add(seq_seq)
                    else:
                        domain_seq_seq[domain] = set()
                        domain_seq_seq[domain].add(seq_seq)

    final_result = open(result_address + 'newpvalue-' + source[8:], 'w')

    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    # file1.readline()
    for line in file1:
        lineSplit = line.rstrip().split("\t")
        dom1 = lineSplit[0]
        dom2 = lineSplit[1]
        AssocScore = lineSplit[10]
        source_score = 0
        if source[8:] == 'intact':
            source_score = lineSplit[2]
        elif source[8:] == 'dip':
            source_score = lineSplit[3]
        elif source[8:] == 'mint':
            source_score = lineSplit[4]
        elif source[8:] == 'biogrid':
            source_score = lineSplit[5]
        elif source[8:] == 'string-exp':
            source_score = lineSplit[6]
        elif source[8:] == 'string-rest':
            source_score = lineSplit[7]
        elif source[8:] == 'sifts_accession':
            source_score = lineSplit[8]
        elif source[8:] == 'hprd':
            source_score = lineSplit[9]

        Ne = 0
        Md = 0
        Kde = 0
        if dom1 in domain_seq_seq:
            Ne = len(domain_seq_seq[dom1])
        if dom2 in domain_seq_seq:
            Md = len(domain_seq_seq[dom2])
        if dom1 in domain_seq_seq and dom2 in domain_seq_seq:
            temp = domain_seq_seq[dom1] & domain_seq_seq[dom2]
            Kde = len(temp)

        if Kde == 0 or source_score == '0':
            final_result.write(str(dom1) + "\t" + str(dom2) + "\t" + str(AssocScore) + "\t" + "NA" + "\n")
            # print("-")
            continue
        Nlist = []
        Nelist = []
        Mdlist = []
        #         Kdelist = []

        NMinusMdlist = []
        NMinusNelist = []

        minMdNe = min(Ne, Md)

        coFactorNe = sqrt(2 * (math.pi) * Ne)
        coFactorMd = sqrt(2 * (math.pi) * Md)
        coFactorN = sqrt(2 * (math.pi) * N)
        NminusNe = N - Ne
        NminusMd = N - Md
        coFactorNminusNe = sqrt(2 * (math.pi) * NminusNe)
        coFactorNminusNMd = sqrt(2 * (math.pi) * NminusMd)
        headCoFactors = coFactorNe * coFactorMd * coFactorNminusNe * coFactorNminusNMd
        headCoFactors = math.log10(headCoFactors)

        logNe = math.log10(Ne / (math.e))
        logNe = logNe * Ne
        logMd = math.log10(Md / (math.e))
        logMd = logMd * Md
        logNMinusNe = math.log10(NminusNe / (math.e))
        logNMinusNe = logNMinusNe * NminusNe
        logNMinusMd = math.log10(NminusMd / (math.e))
        logNMinusMd = logNMinusMd * NminusMd
        headLog = logNe + logMd + logNMinusNe + logNMinusMd

        p_value = 0.0
        for i in range(Kde, minMdNe + 1):
            coFactorI = sqrt(2 * (math.pi) * i)
            NeMinusI = Ne - i
            MdMinusI = Md - i
            NMinusMdNePlusI = N - Md - Ne + i
            if NeMinusI == 0:
                coFactorNeMinusI = 1
            else:
                coFactorNeMinusI = sqrt(2 * (math.pi) * NeMinusI)

            if MdMinusI == 0:
                coFactorMdMinusI = 1
            else:
                coFactorMdMinusI = sqrt(2 * (math.pi) * MdMinusI)
            coFactorNMinusMdNePlusI = sqrt(2 * (math.pi) * NMinusMdNePlusI)

            tailCoFactor = coFactorNeMinusI * coFactorI * coFactorMdMinusI * coFactorNMinusMdNePlusI * coFactorN
            tailCoFactor = math.log10(tailCoFactor)

            logI = math.log10(i / (math.e))
            logI = logI * i
            if NeMinusI == 0:
                logNeMinusI = 1
            else:
                logNeMinusI = math.log10(NeMinusI / (math.e))
                logNeMinusI = logNeMinusI * NeMinusI
            if MdMinusI == 0:
                logMdMinusI = 1
            else:
                logMdMinusI = math.log10(MdMinusI / (math.e))
                logMdMinusI = logMdMinusI * MdMinusI
            logN = math.log10(N / (math.e))
            logN = logN * N
            logNMinusMdNePlusI = math.log10(NMinusMdNePlusI / (math.e))
            logNMinusMdNePlusI = logNMinusMdNePlusI * NMinusMdNePlusI

            tailLog = logI + logNeMinusI + logMdMinusI + logN + logNMinusMdNePlusI
            result = headLog + headCoFactors - tailLog - tailCoFactor
            result = 10 ** result
            p_value += result

        # print(p_value)
        if p_value > 1:
            p_value = "1*"
        final_result.write(str(dom1) + "\t" + str(dom2) + "\t" + str(AssocScore) + "\t" + str(p_value) + "\n")
        continue

    final_result.close()
    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")


def accumulate_pvalues():
    intact_pvalue = dict()
    dip_pvalue = dict()
    mint_pvalue = dict()
    biogrid_pvalue = dict()
    string_exp_pvalue = dict()
    string_rest_pvalue = dict()
    sifts_pvalue = dict()
    hprd_pvalue = dict()
    caps_score = dict()

    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    for line in file1:
        line_sp = line.rstrip().split("\t")
        caps_score[(line_sp[0], line_sp[1])] = line_sp[10]

    file1 = open(result_address + 'newpvalue-intact', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        intact_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]
        # caps_score[(line_sp[0], line_sp[1])] = line_sp[2]

    file1 = open(result_address + 'newpvalue-dip', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dip_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    file1 = open(result_address + 'newpvalue-mint', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        mint_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    file1 = open(result_address + 'newpvalue-biogrid', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        biogrid_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    file1 = open(result_address + 'newpvalue-string-exp', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        string_exp_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    file1 = open(result_address + 'newpvalue-string-rest', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        string_rest_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    file1 = open(result_address + 'newpvalue-sifts', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        sifts_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    file1 = open(result_address + 'newpvalue-hprd', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        hprd_pvalue[(line_sp[0], line_sp[1])] = line_sp[3]

    result = open(result_address + 'newpvalue-all', 'w')
    for item in intact_pvalue:
        result.write(
            item[0] + '\t' + item[1] + '\t' + caps_score[item] + '\t' + intact_pvalue[item] + '\t' + dip_pvalue[item] +
            '\t' + mint_pvalue[item] + '\t' + hprd_pvalue[item] + '\t' + biogrid_pvalue[item] + '\t' +
            string_exp_pvalue[item] + '\t' + string_rest_pvalue[item] + '\t' + sifts_pvalue[item] + '\n')

    result.close()


def gold_silver_bronze():
    gs = set()
    gs_file = open(result_address + 'pfam-pfam-interaction-goldstandard', 'r')
    for line in gs_file:
        line_sp = line.rstrip().split("\t")
        gs.add((line_sp[0], line_sp[1]))
        gs.add((line_sp[1], line_sp[0]))

    calculated_dict = {}
    lenght = 0
    calculated = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    for line in calculated:
        lenght += 1
        line_sp = line.rstrip().split("\t")
        calculated_dict[(line_sp[0], line_sp[1])] = dict()
        calculated_dict[(line_sp[0], line_sp[1])]['caps'] = line_sp[10]
        calculated_dict[(line_sp[0], line_sp[1])]['intact'] = line_sp[2]
        calculated_dict[(line_sp[0], line_sp[1])]['dip'] = line_sp[3]
        calculated_dict[(line_sp[0], line_sp[1])]['mint'] = line_sp[4]
        calculated_dict[(line_sp[0], line_sp[1])]['biogrid'] = line_sp[5]
        calculated_dict[(line_sp[0], line_sp[1])]['string_exp'] = line_sp[6]
        calculated_dict[(line_sp[0], line_sp[1])]['string_rest'] = line_sp[7]
        calculated_dict[(line_sp[0], line_sp[1])]['sifts_acc'] = line_sp[8]
        calculated_dict[(line_sp[0], line_sp[1])]['hprd'] = line_sp[9]

    result = open(result_address + 'result-all', 'w')
    result.write(
        "D1\tD2\tSCORE\tINTACT_SCORE\tINTACT_PV\tDIP_SCORE\tDIP_PV\tMINT_SCORE\tMINT_PV\tHPRD_SCORE\tHPRD_PV"
        "\tBIOGRID_SCORE\tBIOGRID_PV\tSTRING_EXP_SCORE\tSTRING_EXP_PV\tSTRING_REST_SCORE\tSTRING_REST_PV\tSIFTS_SCORE"
        "\tSIFTS_PV\tCLASS\tINTERPRO\n")
    pv = open(result_address + 'newpvalue-all', 'r')
    d_gold_distinct = set()
    d_silver_distinct = set()
    d_bronze_distinct = set()
    for line in pv:
        line_sp = line.rstrip().split("\t")
        d1 = line_sp[0]
        d2 = line_sp[1]
        gs_flag = 'No'
        if (d1, d2) in gs:
            gs_flag = 'Yes'
        caps_score = line_sp[2]
        if caps_score != calculated_dict[(d1, d2)]['caps']:
            print("Something Wrong!")
            print(d1, d2, caps_score, calculated_dict[(d1, d2)]['caps'])
            exit()
        intact = line_sp[3]
        dip = line_sp[4]
        mint = line_sp[5]
        hprd = line_sp[6]
        biogrid = line_sp[7]
        string_exp = line_sp[8]
        string_rest = line_sp[9]
        sifts = line_sp[10]
        critical_val = 0.05 / lenght

        how_many_significant = 0
        how_many = 0

        for i in [intact, dip, mint, hprd, biogrid, string_exp, string_rest, sifts]:
            if i != 'NA':
                how_many += 1
                if i != '1*' and float(i) <= critical_val:
                    how_many_significant += 1

        if how_many >= 4 and how_many == how_many_significant:
            quality = 'Gold'
            d_gold_distinct.add(d1)
            d_gold_distinct.add(d2)
        elif how_many < 4 and how_many == how_many_significant:
            quality = 'Silver'
            d_silver_distinct.add(d1)
            d_silver_distinct.add(d2)
        else:
            quality = 'Bronze'
            d_bronze_distinct.add(d1)
            d_bronze_distinct.add(d2)

        result.write(f"{d1}\t{d2}\t{caps_score}\t{calculated_dict[(d1,d2)]['intact']}\t{intact}\t"
                     f"{calculated_dict[(d1, d2)]['dip']}\t{dip}\t{calculated_dict[(d1,d2)]['mint']}\t{mint}\t"
                     f"{calculated_dict[(d1, d2)]['hprd']}\t{hprd}\t{calculated_dict[(d1, d2)]['biogrid']}\t{biogrid}\t"
                     f"{calculated_dict[(d1, d2)]['string_exp']}\t{string_exp}\t{calculated_dict[(d1,d2)]['string_rest']}"
                     f"\t{string_rest}\t{calculated_dict[(d1, d2)]['sifts_acc']}\t{sifts}\t{quality}\t{gs_flag}\n")

    result.close()

    print('distinct gold doms: ', d_gold_distinct)
    print('distinct silver doms: ', d_silver_distinct)
    print('distinct bronze doms: ', d_bronze_distinct)


def one_to_one():
    result = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'r')
    doms = []
    interactions_dict = dict()
    result.readline()
    for line in result:
        line_sp = line.rstrip().split("\t")
        doms.append(line_sp[0])
        doms.append(line_sp[1])
        if line_sp[0] in interactions_dict:
            interactions_dict[line_sp[0]].add(line_sp[1])
        else:
            interactions_dict[line_sp[0]] = set()
            interactions_dict[line_sp[0]].add(line_sp[1])

        if line_sp[1] in interactions_dict:
            interactions_dict[line_sp[1]].add(line_sp[0])
        else:
            interactions_dict[line_sp[1]] = set()
            interactions_dict[line_sp[1]].add(line_sp[0])

    interactions = set()
    result = open(result_address + 'pfam-pfam-interaction-calculated_tuple', 'r')
    one_to_one = open(result_address + 'one_to_one_tuple', 'w')
    # line = result.readline()
    # one_to_one.write(line)
    for line in result:
        line_sp = line.rstrip().split("\t")
        if len(interactions_dict[line_sp[0]]) == 1 and len(interactions_dict[line_sp[1]]) == 1:
            interactions.add((line_sp[0], line_sp[1]))
            one_to_one.write(line)

    # print((ok_doms))
    print(len(interactions))
