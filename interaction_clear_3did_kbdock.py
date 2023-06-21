import datetime
from main import result_address, source_address


def clean_3did_kbdock_domine_downloaded_files():
    start = datetime.datetime.now()
    print("After Downloading Newest version of 3did and KBDOCK, and DOMAINE data are to be cleaned.")
    file1 = open(source_address + '3did_flat', 'r')
    result = open(result_address + '3did', 'w')
    associations = set()
    for line in file1:
        if 'PF' not in line:
            continue
        line_sp = line.split('\t')
        PF1 = line_sp[len(line_sp) - 2].lstrip()
        PF1 = PF1[1:8]
        PF2 = line_sp[len(line_sp) - 1]
        PF2 = PF2[:7]
        for line in file1:
            if line.startswith('//'):
                break
            if not line.startswith('#=3D'):
                continue
            line_sp = line.split('\t')
            temp = line_sp[2].split(':')
            chain1 = temp[0]
            temp = line_sp[3].split(':')
            chain2 = temp[0]
            if chain1 != chain2:
                if PF1 < PF2:
                    associations.add((PF1, PF2))
                else:
                    associations.add((PF2, PF1))

    for (pf1, pf2) in associations:
        result.write(pf1 + '\t' + pf2 + '\n')

    result.close()

    file1 = open(source_address + 'INTERACTION.txt', 'r')
    result = open(result_address + 'domine', 'w')
    associations = set()
    for line in file1:
        if 'PF' not in line:
            continue
        line_sp = line.split('|')
        PF1 = line_sp[0]
        PF2 = line_sp[1]
        if PF1 < PF2:
            associations.add((PF1, PF2))
        else:
            associations.add((PF2, PF1))

    for (pf1, pf2) in associations:
        result.write(pf1 + '\t' + pf2 + '\n')
    result.close()

    end = datetime.datetime.now()
    print("Running Time: " + str(end - start) + "\n")


def intersection_CAPS_3did_kbdock():
    interactions_3did = set()
    pfam_3did = set()
    file1 = open(result_address + '3did', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        interactions_3did.add((item1, item2))
        pfam_3did.add(item1)
        pfam_3did.add(item2)

    interactions_kbdock = set()
    pfam_kbdock = set()
    # file1 = open(result_address + 'kbdock', 'r')
    # for line in file1:
    #     line_sp = line.rstrip().split('\t')
    #     item1 = line_sp[0]
    #     item2 = line_sp[1]
    #     interactions_kbdock.add((item1, item2))
    #     pfam_kbdock.add(item1)
    #     pfam_kbdock.add(item2)

    interactions_domine = set()
    pfam_domine = set()
    file1 = open(result_address + 'domine', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        interactions_domine.add((item1, item2))
        pfam_domine.add(item1)
        pfam_domine.add(item2)

    interactions_instruct = set()
    pfam_instruct = set()
    # file1 = open(result_address + 'instruct', 'r')
    # for line in file1:
    #     line_sp = line.rstrip().split('\t')
    #     item1 = line_sp[0]
    #     item2 = line_sp[1]
    #     interactions_instruct.add((item1, item2))
    #     pfam_instruct.add(item1)
    #     pfam_instruct.add(item2)

    caps_score = dict()
    caps = set()
    pfam_caps = set()
    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            caps.add((line_sp[0], line_sp[1]))
            caps_score[(line_sp[0], line_sp[1])] = line_sp[9]
        else:
            caps.add((line_sp[1], line_sp[0]))
            caps_score[(line_sp[1], line_sp[0])] = line_sp[9]

        pfam_caps.add(line_sp[0])
        pfam_caps.add(line_sp[1])

    intersecion_idm_3did = caps & interactions_3did
    intersecion_idm_instruct = caps & interactions_instruct
    intersecion_idm_kbdock = caps & interactions_kbdock
    intersecion_idm_domine = caps & interactions_domine
    # for item in intersecion:
    #     print (item)
    #     raw_input()
    print('Instruct:\t' + str(len(interactions_instruct)))
    print('domine:\t' + str(len(interactions_domine)))
    print('IDM & Instruct:\t' + str(len(intersecion_idm_instruct)))
    print('IDM & DOMINE:\t' + str(len(intersecion_idm_domine)))
    print('IDM & kbdock:\t' + str(len(intersecion_idm_kbdock)))
    print('IDM & 3did:\t' + str(len(intersecion_idm_3did)))

    intact_score = dict()
    intact = set()
    pfam_intact = set()
    file1 = open(result_address + 'pfam-pfam-interaction-intact', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            intact.add((line_sp[0], line_sp[1]))
            intact_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            intact.add((line_sp[1], line_sp[0]))
            intact_score[(line_sp[1], line_sp[0])] = line_sp[5]

        pfam_intact.add(line_sp[0])
        pfam_intact.add(line_sp[1])

    dip_score = dict()
    dip = set()
    pfam_dip = set()
    file1 = open(result_address + 'pfam-pfam-interaction-dip', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            dip.add((line_sp[0], line_sp[1]))
            dip_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            dip.add((line_sp[1], line_sp[0]))
            dip_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_dip.add(line_sp[0])
        pfam_dip.add(line_sp[1])

    mint_score = dict()
    mint = set()
    pfam_mint = set()
    file1 = open(result_address + 'pfam-pfam-interaction-mint', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            mint.add((line_sp[0], line_sp[1]))
            mint_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            mint.add((line_sp[1], line_sp[0]))
            mint_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_mint.add(line_sp[0])
        pfam_mint.add(line_sp[1])

    biogrid_score = dict()
    biogrid = set()
    pfam_biogrid = set()
    file1 = open(result_address + 'pfam-pfam-interaction-biogrid', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            biogrid.add((line_sp[0], line_sp[1]))
            biogrid_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            biogrid.add((line_sp[1], line_sp[0]))
            biogrid_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_biogrid.add(line_sp[0])
        pfam_biogrid.add(line_sp[1])

    stringg_exp_score = dict()
    stringg_exp = set()
    pfam_stringg_exp = set()
    file1 = open(result_address + 'pfam-pfam-interaction-string-exp', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            stringg_exp.add((line_sp[0], line_sp[1]))
            stringg_exp_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            stringg_exp.add((line_sp[1], line_sp[0]))
            stringg_exp_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_stringg_exp.add(line_sp[0])
        pfam_stringg_exp.add(line_sp[1])

    stringg_rest_score = dict()
    stringg_rest = set()
    pfam_stringg_rest = set()
    file1 = open(result_address + 'pfam-pfam-interaction-string-rest', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            stringg_rest.add((line_sp[0], line_sp[1]))
            stringg_rest_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            stringg_rest.add((line_sp[1], line_sp[0]))
            stringg_rest_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_stringg_rest.add(line_sp[0])
        pfam_stringg_rest.add(line_sp[1])

    hprd_score = dict()
    hprd = set()
    pfam_hprd = set()
    file1 = open(result_address + 'pfam-pfam-interaction-hprd', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            hprd.add((line_sp[0], line_sp[1]))
            hprd_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            hprd.add((line_sp[1], line_sp[0]))
            hprd_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_hprd.add(line_sp[0])
        pfam_hprd.add(line_sp[1])

    sifts_accession_score = dict()
    sifts_accession = set()
    pfam_sifts_accession = set()
    file1 = open(result_address + 'pfam-pfam-interaction-sifts', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        # if line_sp[0] == 'PF00959' and line_sp[1] == 'PF04965':
        #     print(line_sp)
        #     raw_input()

        if line_sp[1] > line_sp[0]:
            sifts_accession.add((line_sp[0], line_sp[1]))
            sifts_accession_score[(line_sp[0], line_sp[1])] = line_sp[5]
        else:
            sifts_accession.add((line_sp[1], line_sp[0]))
            sifts_accession_score[(line_sp[1], line_sp[0])] = line_sp[5]
        pfam_sifts_accession.add(line_sp[0])
        pfam_sifts_accession.add(line_sp[1])

    sifts = set()
    pfam_sifts = set()
    file1 = open(result_address + 'pfam-pfam-interaction-sifts', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            sifts.add((line_sp[0], line_sp[1]))
        else:
            sifts.add((line_sp[1], line_sp[0]))
        pfam_sifts.add(line_sp[0])
        pfam_sifts.add(line_sp[1])

    print('Database\tPfam-Pfam Interaction')
    print('3did:\t' + str(len(interactions_3did)))
    print('kbdock:\t' + str(len(interactions_kbdock)))
    print('domine:\t' + str(len(interactions_domine)))
    # print('intact & dip: ' + str(len(interactions_caps)))
    print(' ')
    print('IDM result:\t' + str(len(caps)))
    print(' ')
    print('IDM Inputs')
    print('intact:\t' + str(len(intact)))
    print('dip:\t' + str(len(dip)))
    print('mint:\t' + str(len(mint)))
    print('biogrid:\t' + str(len(biogrid)))
    print('string_exp:\t' + str(len(stringg_exp)))
    print('string_rest:\t' + str(len(stringg_rest)))
    print('hprd:\t' + str(len(hprd)))
    print('sifts_accession:\t' + str(len(sifts_accession)))
    print('sifts:\t' + str(len(sifts)))
    print('sifts - sifts_accession:\t' + str(len(sifts - sifts_accession)))
    print(' ')

    intersecion1 = interactions_kbdock & interactions_3did
    print('Observed -> 3did & kbdock:\t' + str(len(intersecion1)))
    gs = (intact | dip | mint | biogrid | stringg_exp | stringg_rest | hprd | sifts_accession) & intersecion1
    print('GS -> Union(IDM-INPUTS) & Observed:\t' + str(len(gs)))
    print(' ')

    intersecion = gs & caps
    print('IDM & GS:\t' + str(len(intersecion)))

    intersecion = gs & interactions_domine
    print('GS & DOMINE:\t' + str(len(intersecion)))

    intersecion = gs & interactions_domine & caps
    print('GS & DOMINE & caps:\t' + str(len(intersecion)))

    intersecion = gs & interactions_instruct
    print('GS & INstruct:\t' + str(len(intersecion)))

    intersecion = interactions_kbdock & caps
    print('IDM & kbdock:\t' + str(len(intersecion)))

    intersecion = interactions_3did & caps
    print('IDM & 3did:\t' + str(len(intersecion)))

    intersecion = gs & intact
    print('GS & intact:\t' + str(len(intersecion)))

    intersecion = gs & dip
    print('GS & dip:\t' + str(len(intersecion)))

    intersecion = gs & mint
    print('GS & mint:\t' + str(len(intersecion)))

    intersecion = gs & biogrid
    print('GS & biogrid:\t' + str(len(intersecion)))

    intersecion = gs & stringg_exp
    print('GS & string_exp:\t' + str(len(intersecion)))

    intersecion = gs & stringg_rest
    print('GS & string_rest:\t' + str(len(intersecion)))

    intersecion = gs & hprd
    print('GS & hprd:\t' + str(len(intersecion)))

    intersecion = gs & sifts_accession
    print('GS & sifts_accession:\t' + str(len(intersecion)))

    print(' ')

    intersecion = caps & intact
    print('IDM & intact:\t' + str(len(intersecion)))

    intersecion = caps & dip
    print('IDM & dip:\t' + str(len(intersecion)))

    intersecion = caps & mint
    print('IDM & mint:\t' + str(len(intersecion)))

    intersecion = caps & biogrid
    print('IDM & biogrid:\t' + str(len(intersecion)))

    intersecion = caps & stringg_exp
    print('IDM & string_exp:\t' + str(len(intersecion)))

    intersecion = caps & stringg_rest
    print('IDM & string_rest:\t' + str(len(intersecion)))

    intersecion = caps & hprd
    print('IDM & hprd:\t' + str(len(intersecion)))

    intersecion = caps & sifts_accession
    print('IDM & sifts_accession:\t' + str(len(intersecion)))

    intersecion = interactions_3did & interactions_kbdock & interactions_domine
    print('3did & kbdock & domine:\t' + str(len(intersecion)))

    intersection = intact & dip
    print('intact & dip:\t' + str(len(intersection)))

    intersection = intact & mint
    print('intact & mint:\t' + str(len(intersection)))

    intersection = intact & biogrid
    print('intact & biogrid:\t' + str(len(intersection)))

    intersection = intact & stringg_exp
    print('intact & string_exp:\t' + str(len(intersection)))

    intersection = intact & stringg_rest
    print('intact & string_rest:\t' + str(len(intersection)))

    intersection = intact & sifts_accession
    print('intact & sifts_accession:\t' + str(len(intersection)))

    intersection = intact & sifts
    print('intact & sifts:\t' + str(len(intersection)))

    intersection = dip & mint
    print('dip & mint:\t' + str(len(intersection)))

    intersection = dip & biogrid
    print('dip & biogrid:\t' + str(len(intersection)))

    intersection = dip & stringg_exp
    print('dip & string_exp:\t' + str(len(intersection)))

    intersection = dip & stringg_rest
    print('dip & string_rest:\t' + str(len(intersection)))

    intersection = dip & sifts_accession
    print('dip & sifts_accession:\t' + str(len(intersection)))

    intersection = dip & sifts
    print('dip & sifts:\t' + str(len(intersection)))

    intersection = mint & biogrid
    print('mint & biogrid:\t' + str(len(intersection)))

    intersection = mint & stringg_exp
    print('mint & string_exp:\t' + str(len(intersection)))

    intersection = mint & stringg_rest
    print('mint & string_rest:\t' + str(len(intersection)))

    intersection = mint & sifts_accession
    print('mint & sifts_accession:\t' + str(len(intersection)))

    intersection = mint & sifts
    print('mint & sifts:\t' + str(len(intersection)))

    intersection = biogrid & stringg_exp
    print('biogrid & string_exp:\t' + str(len(intersection)))

    intersection = biogrid & stringg_rest
    print('biogrid & string_rest:\t' + str(len(intersection)))

    intersection = biogrid & sifts_accession
    print('biogrid & sifts_accession:\t' + str(len(intersection)))

    intersection = biogrid & sifts
    print('biogrid & sifts:\t' + str(len(intersection)))

    intersection = stringg_rest & stringg_exp
    print('string_exp & string_rest:\t' + str(len(intersection)))

    intersection = sifts_accession & stringg_exp
    print('sifts_accession & string_exp:\t' + str(len(intersection)))

    intersection = sifts_accession & stringg_rest
    print('sifts_accession & string_rest:\t' + str(len(intersection)))

    intersection = sifts & stringg_exp
    print('sifts & string_exp:\t' + str(len(intersection)))

    intersection = sifts & stringg_rest
    print('sifts & string_rest:\t' + str(len(intersection)))
    #
    # intersection = sifts & sifts_accession
    # print('sifts & sifts_accession:\t' + str(len(intersection)))
    # print(sifts - sifts_accession)
    #
    # intersection = intact & intersecion1
    # # for line in intersection:
    # #     print(intact_score[line])
    # print('intact & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = dip & intersecion1
    # # for line in intersection:
    # #     print(dip_score[line])
    # print('dip & intersecion1:\t' + str(len(intersection)))
    #
    # intersection = mint & intersecion1
    # # for line in intersection:
    # #     print(mint_score[line])
    # print('mint & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = biogrid & intersecion1
    # # for line in intersection:
    # #     print(biogrid_score[line])
    # print('biogrid & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = stringg_exp & intersecion1
    # print('string_exp & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = stringg_rest & intersecion1
    # print('string_rest & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = sifts_accession & intersecion1
    # # for line in intersection:
    # #     print(sifts_accession_score[line])
    # print('sifts_accession & gold-standard:\t' + str(len(intersection)))
    # # exit()
    #
    # intersection = sifts & intersecion1
    # print('sifts & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = intact & dip & mint & biogrid & stringg_exp & stringg_rest & sifts_accession & intersecion1
    # print('Intersection(CAPS inputs) & gold-standard:\t' + str(len(intersection)))
    #
    # intersection = (intact | dip | mint | biogrid | stringg_exp | stringg_rest | sifts_accession) & intersecion1
    # print('Union(CAPS inputs) & gold-standard:\t' + str(len(intersection)))

    # intersection = sifts_accession & interactions_kbdock
    # print('sifts_accession & interactions_kbdock:\t' + str(len(intersection)))
    #
    # intersection = sifts_accession & interactions_3did
    # print('sifts_accession & interactions_3did:\t' + str(len(intersection)))
    #
    # intersection = sifts & interactions_kbdock
    # print('sifts & interactions_kbdock:\t' + str(len(intersection)))
    #
    # intersection = sifts & interactions_3did
    # print('sifts & interactions_3did:\t' + str(len(intersection)))

    print(' ')
    print(' ')

    print('Distinct Pfams')
    print('pfam kbdock:\t' + str(len(pfam_kbdock)))
    print('pfam 3did:\t' + str(len(pfam_3did)))
    print('pfam domine:\t' + str(len(pfam_domine)))
    print('pfam IDM:\t' + str(len(pfam_caps)))
    print('pfam intact:\t' + str(len(pfam_intact)))
    print('pfam dip:\t' + str(len(pfam_dip)))
    print('pfam mint:\t' + str(len(pfam_mint)))
    print('pfam biogrid:\t' + str(len(pfam_biogrid)))
    print('pfam string_exp:\t' + str(len(pfam_stringg_exp)))
    print('pfam string_rest:\t' + str(len(pfam_stringg_rest)))
    print('pfam sifts_accession:\t' + str(len(pfam_sifts_accession)))
    print('pfam pfam_sifts:\t' + str(len(pfam_sifts)))

    intersection1 = pfam_3did & pfam_kbdock
    print('Observed -> pfam_3did & pfam_kbdock:\t' + str(len(intersection1)))

    #
    # intersection = pfam_intact & pfam_kbdock
    # print('pfam_intact & pfam_kbdock:\t' + str(len(intersection)))
    #
    # intersection = pfam_dip & pfam_kbdock
    # print('pfam_dip & pfam_kbdock:\t' + str(len(intersection)))
    #
    # intersection = pfam_mint & pfam_kbdock
    # print('pfam_mint & pfam_kbdock:\t' + str(len(intersection)))
    #
    # intersection = pfam_biogrid & pfam_kbdock
    # print('pfam_biogrid & pfam_kbdock:\t' + str(len(intersection)))
    #
    # intersection = pfam_stringg_exp & pfam_kbdock
    # print('pfam_stringg & pfam_kbdock:\t' + str(len(intersection)))
    #
    # intersection = pfam_sifts_accession & pfam_kbdock
    # print('pfam_sifts_accession & pfam_kbdock:\t' + str(len(intersection)))
    #
    # intersection = pfam_sifts & pfam_kbdock
    # print('pfam_sifts & pfam_kbdock:\t' + str(len(intersection)))
    # dif = pfam_kbdock - pfam_sifts
    # print(dif)
    #
    # intersection = pfam_intact & pfam_3did
    # print('pfam_intact & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_dip & pfam_3did
    # print('pfam_dip & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_mint & pfam_3did
    # print('pfam_mint & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_biogrid & pfam_3did
    # print('pfam_biogrid & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_stringg_exp & pfam_3did
    # print('pfam_stringg & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_sifts_accession & pfam_3did
    # print('pfam_sifts_accession & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_sifts & pfam_3did
    # print('pfam_sifts & pfam_3did:\t' + str(len(intersection)))
    #
    # intersection = pfam_intact & intersection1
    # print('pfam_intact & GS:\t' + str(len(intersection)))
    #
    # intersection = pfam_dip & intersection1
    # print('pfam_dip & GS:\t' + str(len(intersection)))
    #
    # intersection = pfam_mint & intersection1
    # print('pfam_mint & GS:\t' + str(len(intersection)))
    #
    # intersection = pfam_biogrid & intersection1
    # print('pfam_biogrid & GS:\t' + str(len(intersection)))
    #
    # intersection = pfam_stringg_exp & intersection1
    # print('pfam_stringg & GS:\t' + str(len(intersection)))
    #
    # intersection = pfam_sifts_accession & intersection1
    # print('pfam_sifts_accession & GS:\t' + str(len(intersection)))
    #
    # intersection = pfam_sifts & intersection1
    # print('pfam_sifts & GS:\t' + str(len(intersection)))

    intersection = pfam_sifts & pfam_sifts_accession & pfam_stringg_exp & pfam_biogrid & pfam_mint & pfam_dip & pfam_intact & intersection1
    print('Intersection(IDM inputs) & Observed:\t' + str(len(intersection)))

    gs = (
                 pfam_sifts | pfam_sifts_accession | pfam_stringg_exp | pfam_biogrid | pfam_mint | pfam_dip | pfam_intact) & intersection1
    print('GS -> Union(IDM inputs) & Observed:\t' + str(len(gs)))

    intersection = gs & pfam_caps
    print('caps & GS:\t' + str(len(intersection)))

    # for item in intersecion3:
    #     item = item[0]+item[1]
    #     if item in caps_score:
    #         print(caps_score[item])


def kbdock_union_3did():
    interactions_3did = set()
    pfam_3did = set()
    file1 = open(result_address + '3did', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        interactions_3did.add((item1, item2))
        pfam_3did.add(item1)
        pfam_3did.add(item2)

    interactions_kbdock = set()
    pfam_kbdock = set()
    file1 = open(result_address + 'kbdock', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        interactions_kbdock.add((item1, item2))
        pfam_kbdock.add(item1)
        pfam_kbdock.add(item2)

    union = interactions_kbdock | interactions_3did
    print(len(union))
    result = open(result_address + 'unionkbdock3did', 'w')
    for item in union:
        flag = ''
        if item in interactions_kbdock and item in interactions_3did:
            flag = 'both'
        elif item in interactions_kbdock:
            flag = 'kbdock'
        elif item in interactions_3did:
            flag = '3did'
        result.write(item[0] + '\t' + item[1] + '\t' + flag + '\n')

    result.close()


def checkOneOnedomain():
    did = set()
    file1 = open(result_address + '3did', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            did.add((item1, item2))
        else:
            did.add((item2, item1))

    kbdock = set()
    file1 = open(result_address + 'kbdock', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        item1 = line_sp[0]
        item2 = line_sp[1]
        if item1 < item2:
            kbdock.add((item1, item2))
        else:
            kbdock.add((item2, item1))

    intersection = kbdock & did
    print(len(did))
    print(len(kbdock))
    print(len(intersection))

    caps_score = dict()
    caps = set()
    file1 = open(result_address + 'pfam-pfam-interaction-calculated', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if line_sp[1] > line_sp[0]:
            caps.add((line_sp[0], line_sp[1]))
            caps_score[(line_sp[0], line_sp[1])] = line_sp[9]
        else:
            caps.add((line_sp[1], line_sp[0]))
            caps_score[(line_sp[1], line_sp[0])] = line_sp[9]

    print(len(caps))

    all_interactions = set()
    for source in ['source1_intact', 'source2_mint', 'source3_dip', 'source4_biogrid', 'source5_string-exp',
                   'source5_string-rest', 'source6_sifts', 'source7_hprd']:
        # for source in ['source5_string-exp', 'source5_string-rest']:
        file1 = open(result_address + source, 'r')
        for line in file1:
            if '_' in line:
                continue
            line_sp = line.rstrip().split('\t')
            seq1 = line_sp[0]
            seq2 = line_sp[1]
            all_interactions.add((seq1, seq2))
    print(len(all_interactions))

    seqDom = dict()
    file1 = open(result_address + 'pfam-seq-sp', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[0]
        seq = line_sp[1]
        if seq in seqDom:
            seqDom[seq].append(dom)
        else:
            seqDom[seq] = list()
            seqDom[seq].append(dom)

    file1 = open(result_address + 'pfam-seq-tr', 'r')
    for line in file1:
        line_sp = line.rstrip().split('\t')
        dom = line_sp[0]
        seq = line_sp[1]
        if seq in seqDom:
            seqDom[seq].append(dom)
        else:
            seqDom[seq] = list()
            seqDom[seq].append(dom)

    counter = set()
    counter2 = set()
    pf_lst = set()
    seq_lst = set()
    for (seq1, seq2) in all_interactions:
        if seq1 in seqDom and seq2 in seqDom and len(seqDom[seq1]) == 1 and len(seqDom[seq2]) == 1:
            if seq1 < seq2:
                counter.add((seq1, seq2))
            else:
                counter.add((seq2, seq1))

            seq_lst.add(seq1)
            seq_lst.add(seq2)

            pf1 = seqDom[seq1][0]
            pf2 = seqDom[seq2][0]
            if (pf1, pf2) in caps or (pf2, pf1) in caps:
                if pf1 < pf2:
                    counter2.add((pf1, pf2))
                else:
                    counter2.add((pf2, pf1))
                pf_lst.add(pf1)
                pf_lst.add(pf2)
                # print(seq1, seq2, pf1, pf2)
                # raw_input()

    print(len(counter))
    print(len(seq_lst))
    print(len(counter2))
    print(len(pf_lst))
    print(len(counter2 & intersection))

    res_file = open(result_address + 'result', 'w')
    file1 = open(result_address + 'result-all', 'r')
    line = file1.readline()
    res_file.write(line.rstrip() + "\tSINGLE_DOMAIN_SEQUENCES\n")
    for line in file1:
        line_sp = line.rstrip().split('\t')
        if (line_sp[0], line_sp[1]) in counter2 or (line_sp[1], line_sp[0]) in counter2:
            res_file.write(line.rstrip() + "\tSingle\n")
        else:
            res_file.write(line.rstrip() + "\tMultiple\n")

    res_file.close()

    print(len(caps))
