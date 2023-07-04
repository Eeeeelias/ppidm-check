import pickle
import sys
import timeit

result_address = 'resultdata/'
source_address = 'sourcedata/'


import interaction_clear_3did_kbdock as ic3k
import process_tables as pt
import filtering
import pvalue as pv


def similarity_calculator_interaction_wrapper(args):
    source, seqDom, seqpdbchain, pdbchainDom = args
    pt.similarity_calculator_interaction(source, 'pfam', seqDom, seqpdbchain, pdbchainDom)


if __name__ == '__main__':
    sources = ['source1_intact', 'source2_mint', 'source3_dip', 'source4_biogrid', 'source5_string-exp',
               'source5_string-rest', 'source6_sifts', 'source7_hprd']

    start = timeit.default_timer()
    # seqDom, seqpdbchain, pdbchainDom = pt.read_chain_dom()
    # pickle.dump((seqDom, seqpdbchain, pdbchainDom), open('seqDom_seqpdbchain_pdbchainDom.pickle', 'wb'))
    # seqDom, seqpdbchain, pdbchainDom = pickle.load(open('pickles/seqDom_seqpdbchain_pdbchainDom.pickle', 'rb'))
    # print("Loading files from pickle took:", round(timeit.default_timer() - start, 1), "seconds")

    # for i in sources:
    #      pt.similarity_calculator_interaction(i, 'pfam', seqDom, seqpdbchain, pdbchainDom)
    #
    # # sifts_reader_process('sifts', 'pfam')
    # ic3k.clean_3did_kbdock_domine_downloaded_files()
    #
    # filtering.create_wrong_assocations()
    #
    # filtering.assign_interaction()
    #
    # ic3k.intersection_CAPS_3did_kbdock()
    # ic3k.kbdock_union_3did()

    # for i in sources:
    #      pv.pvalue_calculation(i, seqDom, pdbchainDom)
    pv.accumulate_pvalues()
    pv.gold_silver_bronze()

    # pv.one_to_one()
    print("Took:", timeit.default_timer() - start)
    # checkOneOnedomain()
    # verify_go_for_each_pair()

