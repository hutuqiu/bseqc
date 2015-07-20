#!/usr/bin/env python

'''
Description:  BSeQC main executable function

Input:  the standard SAM files
Output: the bias-free SAM files, Mbias Plot/Table, Duplicate read distribution, Filter Report
Support:
1. Use M-bias plot to measure the BS-specific biases (eg: overhang end-repaired bias; 5 prime non-conversion bias)
   and sequence related biases (eg: sequence into adaptor; low-quality nucleotide call).
2. Trim the biases based on the M-bias plot automatically
3. Remove the duplicate reads resulting from possible over-amplification
4. Keep only one copy of the overlapping segment of two read mates in paired-end seq
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
from BSeQC.qc_assess import qc_report as QR
from BSeQC.qc_filter import qc_filter as QF
from BSeQC.read import check_file as check

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

    if len(options.trim_file) != 0 and not os.path.isfile(options.trim_file):
        error("Can't open the ref file: " + options.trim_file)

    options.read_length = options.read_length.split(',')
    sam_inf = options.sam_file
    ref_file = options.ref_file
    bsm = options.bsm
    s_path = options.samtools
    name = options.name
    read_l = options.read_length
    auto = options.automatically
    pvalue = options.pvalue
    drift = options.drift
    trim_file = options.trim_file
    remove_overlap = options.remove_overlap
    filter_dup = options.filter_dup
    p_poisson = options.p_poisson
    gsize = options.gsize
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

    if filter_dup:
        ## if filter_up is TRUE, the duplicate reads will be assessed and shown in Dup_dis.pdf
        ## and the loc_dict & max_cov will be used in the trimming step
        if len(trim_file) != 0:
            info("The trimming file has been defined. But the filter_dup has been set True.")
            info("QC_report will just generate Dup distribution!!")
            info("And the user defined trimming file will be used in the trimming step!!")
        else:
            info("The filter_dup has been set True.")
            info("QC_report not only includes Mbias plot, Mbias table and trimming file, but also Dup distribution.")
        QC_report_MD = QR.QC_Report_Mbias_Dup(sam_inf, ref_file, bsm, s_path, name, read_l, single_on, pvalue, drift,
                                              trim_file,
                                              p_poisson, gsize)
        strand_t, loc_dict, max_cov = QC_report_MD.generator()

    else:
        if len(trim_file) != 0:
            info("The trimming file has been defined. So Ignore the ")
        info("The filter_dup has been set False!! QC_report only includes Mbias plot, Mbias table and trimming file.")
        info("And ignore the collection of the location information for removing duplicate reads!!")
        QC_report_M = QR.QC_Report_Mias(sam_inf, ref_file, bsm, s_path, name, read_l, single_on, pvalue, drift,
                                        trim_file)
        strand_t = QC_report_M.generator()
        #no duplicate location information
        loc_dict = {}
        max_cov = 10000

    if ((auto or filter_dup) and single_on) or ((auto or filter_dup or remove_overlap) and not single_on):
    ## for single-end: qc_filter Mbias or filter duplicate reads
    ## for paired-end: qc_filter Mbias, keep one copy of the overlapping segment, or filter duplicate reads
        info("Start to filter read...")
        if auto:
            info("Automatically trim Mbias...")
        else:
            info("--auto has been set %s ! Ignore trimming Mbias!!" % auto)
        if filter_dup:
            info("Filter duplicate reads...")
        else:
            info("--filter_dup has been set %s ! Ignore removing duplicate reads!!" % filter_dup)
        if remove_overlap and not single_on:
            info("Keep one copy of the overlapping segment...")
        if not remove_overlap and not single_on:
            info(
                "--remove_overlap has been set %s ! Ignore removing one copy of the overlapping segment!!" % remove_overlap)
        if not_mapping:
            info("Keep the not_unique mapping reads!")
        else:
            info("Remove the not_unique mapping reads!!")
        QF.filter_sam(sam_inf, ref_file, bsm, strand_t, read_l, single_on, name, s_path, auto, remove_overlap, loc_dict, max_cov,
                      not_mapping)
        info("Get the filtered SAM file!")
    else:
        if single_on:
            info("Skip the trimming Mbias and removing duplicate reads!!")
            info("Not BSeQC filter report!!")
        else:
            info("Skip the trimming Mbias, removing duplicate reads and removing one copy of the overlapping segment!!")
            info("Not BSeQC filter report!!")

