#!/usr/bin/python


# ------------------------------------
#python package
# ------------------------------------
import sys
import logging

# ------------------------------------
#own python modules
# ------------------------------------
from BSeQC.read.read_info import MappingReader as RI

# ------------------------------------
#logging object
# ------------------------------------

logging.basicConfig(level=20,
                    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w")
info = logging.info
error = logging.error
warning = logging.warning


def get_meth_state(strand, seq, refseq, index):
    '''
    get the methylation state for a specific C
    '''
    BS_conversion = {'+': ('C', 'T'), '-': ('G', 'A')}
    match, convert = BS_conversion[strand[0]]
    ref_index = index + 1
    if strand[1] == '+':
        context = refseq[ref_index: ref_index + 2]
    else:
        context = refseq[ref_index - 1: ref_index + 1]
    if (strand == '++' or strand == '--') and context == 'CG':
        meth_all = 1
    elif (strand == '-+' or strand == '+-') and context == 'GC':
        meth_all = 1
    else:
        meth_all = 0
    if meth_all == 1:
        if seq[index] == convert:
            meth_state = 0
        else:
            meth_state = 1
    else:
        meth_state = 0
    return meth_state, meth_all


def record_site(read, ref, bsm, dige_site):
    '''
    record the location of restriction enzyme digestion site in the read,
    and the methylation level of the restriction enzyme digestion site
    :param read:
    :param ref:
    :param dige_site:
    '''


    #Using the read class
    #If the read is unique mapping or unique and  paired mapping, we will get a information list.
    #If not, we will get a empty list ([]).
    #In: single unique mapping read  Out: [flag,strand,chr,pos,CIGAR,seq,score]
    #In: paired unique mapping read  Out: [flag,strand,chr,pos1,CIGAR,pos2,insert,seq,score]
    read_s = RI(read, bsm)
    read_info = read_s.extract_information()
    strand = ''
    site_meth_list = []
    reverse_strand = ['-+', '+-']

    if len(read_info) == 0:
        return read_info, site_meth_list


    strand, chr, pos, seq = read_info[1], read_info[2], int(read_info[3]) - 1, read_info[-2]
    #print strand, chr, pos, seq

    if len(chr.split('_')) > 1:   #filter some 'chr_' chromosome
        return read_info, site_meth_list

    readlen = len(seq)
    refseq = ref[chr][pos - 1: pos + readlen + 1]


    if strand not in reverse_strand:
    #get the location of restriction enzyme digestion site
        dige_5_index = refseq.find('CCGG')
        dige_3_index = refseq.rfind('CCGG')

    else:
    # if the strand is '-+' (single-end) or '+-', reverse the seq and refseq
        seq = seq[::-1]
        refseq = refseq[::-1]
        dige_5_index = refseq.find('GGCC')
        dige_3_index = refseq.rfind('GGCC')

    if dige_3_index == -1 or dige_5_index == -1:
        warning("Can't find MspI site in the read: %s" %read.rstrip())
        #return dige_dict
        return read_info, site_meth_list


    # get the location of the restriction enzyme digestion site
    #dige_dict[strand]['s'][0].append(dige_5_index + 1)
    #dige_dict[strand]['e'][0].append(dige_3_index + 1)


    # get the methylation state of the restriction enzyme digestion site
    if strand[1] == '-':
        #stand: +- or --, check the third nucleotide of C-CGG
        dige_5_meth, dige_5_all = get_meth_state(strand, seq, refseq, dige_5_index + 1)
        dige_3_meth, dige_3_all = get_meth_state(strand, seq, refseq, dige_3_index + 1)
    else:
        #stand: ++ or -+, check the second nucleotide of C-CGG
        dige_5_meth, dige_5_all = get_meth_state(strand, seq, refseq, dige_5_index)
        dige_3_meth, dige_3_all = get_meth_state(strand, seq, refseq, dige_3_index)



    #dige_dict[strand]['s'][1] = [sum(x) for x in zip(dige_dict[strand]['s'][1], [dige_5_meth, dige_5_all])]
    #dige_5_meth_list = [dige_5_meth] + [0] * 10
    #dige_5_all_list = [dige_5_all] + [0] * 10
    #dige_3_meth_list = [0] * 11
    #dige_3_all_list = [0] * 11
    #calculate the methylation states of 10 nucleotides after dige_5_index
    #for i in range(1, 11):
        #meth_5, all_5 = get_meth_state(strand, seq, refseq, dige_5_index + 2 + i)  # +2 to scan the nucleotide after CCGG
        #meth_3, all_3 = get_meth_state(strand, seq, refseq, len(seq) - i)
        #dige_5_meth_list[i] = meth_5
        #dige_5_all_list[i] = all_5
        #dige_3_meth_list[10 - i] = meth_3
        #dige_3_all_list[10 - i] = all_3
    #site_meth_list = [dige_5_meth_list, dige_5_all_list, dige_3_meth_list, dige_3_all_list, dige_5_index, len(seq)]
    site_meth_list = [dige_5_meth, dige_5_all, 0, 0, dige_5_index, len(seq)]
    if dige_3_index != 0 and dige_3_index != dige_5_index:
        if strand[1] == '+':
            #check the 3' digestion site for '++' and '-+'
            if len(seq[(dige_3_index - 1):(dige_3_index + 3)]) == 4:
                if (strand == '++' and seq[dige_3_index + 2] == 'A') or (strand == '-+' and seq[dige_3_index + 2] == 'T'):
                    site_meth_list = [dige_5_meth, dige_5_all, dige_3_meth, dige_3_all, dige_5_index, dige_3_index]
                    #dige_dict[strand]['e'][1] = [sum(x) for x in zip(dige_dict[strand]['e'][1], [dige_3_meth, dige_3_all])]
                    #dige_dict[strand]['e'] = [sum(x) for x in zip(dige_dict[strand]['e'], [dige_3_meth, dige_3_all])]
                    #for i in range(1, 11):
                    #    meth_3, all_3 = get_meth_state(strand, seq, refseq, dige_3_index - 1 - i)
                    #    dige_3_meth_list[10 - i] = meth_3
                    #    dige_3_all_list[10 - i] = all_3
                    #dige_3_meth_list[-1] = dige_3_meth
                    #dige_3_all_list[-1] = dige_3_all
                    #site_meth_list = [dige_5_meth_list, dige_5_all_list, dige_3_meth_list, dige_3_all_list, dige_5_index, dige_3_index]
                    #print "end_repair_3", refseq[dige_3_index:(dige_3_index + 4)], dige_3_meth_list, dige_3_all_list


        else:
            site_meth_list = [dige_5_meth, dige_5_all, dige_3_meth, dige_3_all, dige_5_index, dige_3_index]
            #dige_dict[strand]['e'][1] = [sum(x) for x in zip(dige_dict[strand]['e'][1], [dige_3_meth, dige_3_all])]
            #dige_dict[strand]['e'] = [sum(x) for x in zip(dige_dict[strand]['e'], [dige_3_meth, dige_3_all])]
            #for i in range(1, 11):
            #    meth_3, all_3 = get_meth_state(strand, seq,refseq, dige_3_index - 1 - i)
            #    dige_3_meth_list[10 - i] = meth_3
            #    dige_3_all_list[10 - i] = all_3
            #dige_3_meth_list[-1] = dige_3_meth
            #dige_3_all_list[-1] = dige_3_all
            #site_meth_list = [dige_5_meth_list, dige_5_all_list, dige_3_meth_list, dige_3_all_list, dige_5_index, dige_3_index]
            #print 'bbbb', dige_3_meth_list, dige_3_all_list

    #return dige_dict
    return read_info, site_meth_list




