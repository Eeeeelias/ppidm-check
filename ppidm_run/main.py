import os
import pickle
import random
import sys
import timeit

result_address = '../resultdata/'
source_address = '../sourcedata/'


import ppidm_run.interaction_clear_3did_kbdock as ic3k
import ppidm_run.process_tables as pt
import ppidm_run.filtering as filtering
import ppidm_run.pvalue as pv


if __name__ == '__main__':
    sources = [x for x in os.listdir(source_address) if x.startswith('source')]
    print("Sources:", sources)

    start = timeit.default_timer()
    # This part gets all the domain-protein information (what proteins are associated with which domains etc.)
    # seqDom, seqpdbchain, pdbchainDom = pt.read_chain_dom()
    # pickle.dump((seqDom, seqpdbchain, pdbchainDom), open('../pickles/seqDom_seqpdbchain_pdbchainDom.pickle', 'wb'))
    seqDom, seqpdbchain, pdbchainDom = pickle.load(open('../pickles/seqDom_seqpdbchain_pdbchainDom.pickle', 'rb'))
    # print("Loading files from pickle took:", round(timeit.default_timer() - start, 1), "seconds")

    # This calculates all the similarity scores for each source
    for i in sources:
         pt.similarity_calculator_interaction(i, 'pfam', seqDom, seqpdbchain, pdbchainDom, redo=False)
    #
    # # # sifts_reader_process('sifts', 'pfam')
    # # This is a cleanup of the domain interaction sources
    ic3k.clean_3did_kbdock_domine_downloaded_files()
    #
    # # This function creates random wrong associations (retaining node degree) for the ddi inference
    filtering.create_wrong_assocations(sources)

    # This function assigns the interactions. It also does all the "hyperparameter optimization"
    filtering.assign_interaction(sources)

    # ic3k.kbdock_union_3did()
    #
    # # This part calculates and sorts all the p-values for every domain-domain interaction
    # for i in sources:
    #      pv.pvalue_calculation(i, seqDom, pdbchainDom)
    # pv.accumulate_pvalues(sources)
    # pv.gold_silver_bronze()

    # pv.one_to_one()
    print("Took:", timeit.default_timer() - start)
    # checkOneOnedomain()
    # verify_go_for_each_pair()

