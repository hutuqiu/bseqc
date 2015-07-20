#!/usr/bin/env python

'''
Description:  Alternative executable function for RRBS in BSeQC

Input:  the standard SAM files
Output: the bias-free SAM files
Support:
1. Use restriction enzyme digestion sites(C-CGG) to remove end-repaired bias
2. Keep only one copy of the overlapping segment of two read mates in paired-end seq

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@author: Xueqiu Lin
@contact: xueqiu.lin@gmail.com
'''

# ------------------------------------
# python modules
# ------------------------------------
import os
import sys
import logging

# ------------------------------------
# own python modules
# ------------------------------------
from BSeQC.read import check_file as check
from BSeQC.read import get_reference as GR
from BSeQC.qc_assess import rrbs_report as RR
from BSeQC.read import dige_information as DI
from BSeQC.qc_filter import rrbs_filter as RF


# ------------------------------------
#logging object
# ------------------------------------

logging.basicConfig(level=20,
                    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
)
info = logging.info
error = logging.error
warning = logging.warning


def run(args):
    options = args.parse_args()
    if len(options.sam_file) == 0:
        error("Missing the SAM file, use -s or --sam option.")
    else:
        options.sam_file = options.sam_file.split(',')
    for s in options.sam_file:
        if not os.path.isfile(s):
            error("Can't open the SAM file: " + s)
            sys.exit(1)

    if len(options.ref_file) == 0:
        error("Missing the reference genome fasta file, use -r or --ref option.")
    else:
        if not os.path.isfile(options.ref_file):
            error("Can't open the ref file: " + options.ref_file)

    if len(options.samtools) != 0:
        if options.samtools[-1] != '/':
            options.samtools += '/'

    if len(options.name) == 0:
        error("Missing the output file name, use -n or --name options.")

    sam_inf = options.sam_file
    ref_file = options.ref_file
    bsm = options.bsm
    s_path = options.samtools
    name = options.name
    dige_site = options.dige_site
    remove_overlap = options.remove_overlap
    not_mapping = options.not_mapping


    info("Get the all parameter!!")
    #check the input mapping files
    sam_format, read_inf = check.check_mapping_file_flag(sam_inf[0], s_path)
    pre_flag = read_inf.readline().split('\t')[1]
    if 'p' in pre_flag:
        single_on = False
        info("The input mapping files are paired-end sequencing!")
    else:
        single_on = True
        info("The input mapping files are single-end sequencing!")

    #get reference information
    ref = GR.get_ref(ref_file)

    ##scan MspI site and trim the end-repaired C
    dige_dict, all_reads, all_mapping_bp, not_mapping_reads, filter_not_mapping_reads, filter_MspI_endrepair_bp, filter_remove_overlap_bp = parser_trim_sambam(
        sam_inf, ref, bsm, s_path, dige_site, single_on, remove_overlap, not_mapping, name)


    ##produce MspI Mbias plot
    RR.generator(dige_dict, single_on, name)

    ##produce the filter report
    report(all_reads, all_mapping_bp, not_mapping_reads, filter_not_mapping_reads, filter_MspI_endrepair_bp,
           filter_remove_overlap_bp, single_on, name)


def parser_trim_sambam(sam_inf, ref, bsmp, s_path, dige_site, single_on, remove_overlap, not_mapping, name):
    '''
    1. Scan each SAM file to calculate the methylation level of the end-repaired C in the MspI site.
    2. Trim the end-repaired C and adapter sequence by recognizing the MspI site
    '''

    info('Scan each SAM file to calculate the methylation level of the end-repaired C in the MspI site.')
    info('and trim the end-repaired C and adapter sequence by recognizing the MspI site.')
    dige_dict = {}
    #dige_dict['++'] = {'s': [[], [0, 0]], 'e': [[], [0, 0]]}
    #dige_dict['-+'] = {'s': [[], [0, 0]], 'e': [[], [0, 0]]}
    #dige_dict['++'] = {'s': [0] * 22, 'e': [0] * 22}
    #dige_dict['-+'] = {'s': [0] * 22, 'e': [0] * 22}
    dige_dict['++'] = {'s': [0] * 2, 'e': [0] * 2}
    dige_dict['-+'] = {'s': [0] * 2, 'e': [0] * 2}
    if not single_on:
        #dige_dict['+-'] = {'s': [[], [0, 0]], 'e': [[], [0, 0]]}
        #dige_dict['--'] = {'s': [[], [0, 0]], 'e': [[], [0, 0]]}
        #dige_dict['+-'] = {'s': [0] * 22, 'e': [0] * 22}
        #dige_dict['--'] = {'s': [0] * 22, 'e': [0] * 22}
        dige_dict['+-'] = {'s': [0] * 2, 'e': [0] * 2}
        dige_dict['--'] = {'s': [0] * 2, 'e': [0] * 2}

    filter_remove_overlap_bp = 0
    filter_not_mapping_reads = 0
    all_reads = 0
    not_mapping_reads = 0
    all_mapping_bp = 0
    filter_MspI_end_repaired_bp = 0

    for s in range(len(sam_inf)):
        sam_format, read_inf = check.check_mapping_file_header(sam_inf[s], s_path)
        out_sam = sam_inf[s][:-4] + '_' + name + '_filter.sam'
        out = open(out_sam, 'w')
        record_mate = {}
        for read in read_inf:
            if read.startswith('@'):
                out.write(read)
                continue
            else:
                all_reads += 1
                #get read information and MspI site
            read_info, site_meth_list = DI.record_site(read, ref, bsmp, dige_site)

            if len(read_info) == 0:
                not_mapping_reads += 1
                if not_mapping:
                    out.write(read)
                else:
                    filter_not_mapping_reads += 1
                continue

            strand = read_info[1]
            if single_on:
                all_mapping_bp += len(read_info[5])
            else:
                all_mapping_bp += len(read_info[7])

            if len(site_meth_list) != 0:
            #record the methylation state of the MspI site
                dige_dict[strand]['s'] = [sum(x) for x in zip(dige_dict[strand]['s'], site_meth_list[:2])]
                dige_dict[strand]['e'] = [sum(x) for x in zip(dige_dict[strand]['e'], site_meth_list[2:])]
                #dige_dict[strand]['s'] = [sum(x) for x in zip(dige_dict[strand]['s'], site_meth_list[0] + site_meth_list[1])]
                #dige_dict[strand]['e'] = [sum(x) for x in zip(dige_dict[strand]['e'], site_meth_list[2] + site_meth_list[3])]

                #trim the end-repaired nucleotides and adapter sequence by the MspI site
                record_mate, filter_MspI_end_repaired_bp, filter_remove_overlap_bp = RF.rrbs_trim(read, read_info,
                                                                                                  site_meth_list,
                                                                                                  single_on,
                                                                                                  remove_overlap,
                                                                                                  record_mate,
                                                                                                  filter_MspI_end_repaired_bp,
                                                                                                  filter_remove_overlap_bp, out)
        out.close()
    return dige_dict, all_reads, all_mapping_bp, not_mapping_reads, filter_not_mapping_reads, filter_MspI_end_repaired_bp, filter_remove_overlap_bp


def report(all_reads, all_mapping_bp, not_mapping_reads, filter_not_mapping_reads, filter_MspI_endrepaire_bp,
           filter_remove_overlap_bp, single_on, name):
    info('Produce the report file...')
    report_out = open(name + "_BSeQC_rrbs_filter_report.txt", 'w')
    report_out.write('Total reads: %d\n' % all_reads)
    if single_on:
        report_out.write('Not unique mapping reads: %d(%.2f%s all reads)\n' % (
            not_mapping_reads, float(not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('Unique mapping reads: %d(%.2f%s all reads)\n' % (
            (all_reads - not_mapping_reads), float(all_reads - not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('Skip not unique mapping reads: %d(%.2f%s all reads)\n' % (
            filter_not_mapping_reads, float(filter_not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('In unique mapping reads:\n')
        report_out.write('All unique mapping basepairs: %d\n' % all_mapping_bp)
        report_out.write("Filter end-repaired nucleotides in MspI site: %d(%.2f%s of unique mapping basepairs)\n" % (
            filter_MspI_endrepaire_bp, float(filter_MspI_endrepaire_bp) / all_mapping_bp * 100, "%"))
    else:
        report_out.write('Not unique paired mapping reads: %d(%.2f%s)\n' % (
            not_mapping_reads, float(not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('Unique paired mapping reads: %d(%.2f%s)\n' % (
            (all_reads - not_mapping_reads), float(all_reads - not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('Skip not paired unique mapping reads: %d(%.2f%s)\n' % (
            filter_not_mapping_reads, float(filter_not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('In unique paired mapping reads:\n')
        report_out.write('All unique paired mapping basepairs: %d\n' % all_mapping_bp)
        report_out.write("Filter end-repaired nucleotides in MspI site : %d(%.2f%s of unique mapping basepairs)\n" % (
            filter_MspI_endrepaire_bp, float(filter_MspI_endrepaire_bp) / all_mapping_bp * 100, "%"))
        report_out.write('Filter overlapped basepairs: %d(%.2f%s of unique paired mapping basepairs)\n' % (
            filter_remove_overlap_bp, float(filter_remove_overlap_bp) / all_mapping_bp * 100, "%"))
    report_out.close()
    info('Get the report file!')



