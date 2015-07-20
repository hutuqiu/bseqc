#!/usr/bin/env python

'''
Description:  Alternative executable function in BSeQC

Input:  the standard SAM files
Output: the bias-free SAM files
Support:
1. Use the strategy in Bis-SNP to trim 5' bisulfite conversion failures
   (For each read, we walk along the read from 5' to 3', and we remove any Cs on the C-strand until we reach the first reference C which is converted to a T)
2. Remove the duplicate reads resulting from possible over-amplification
3. Keep only one copy of the overlapping segment of two read mates in paired-end seq

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
from BSeQC.read.read_info import MappingReader as RI
from BSeQC.qc_filter import nonuniform_filter as NF
from BSeQC.read import check_file as check
from BSeQC.read import Loc_information as LI
from BSeQC.qc_filter.duplicate_filter import duplicate_filter as DF
from BSeQC.qc_assess import duplicate_report as DR
from BSeQC.read import get_reference as GR
from BSeQC.qc_assess import nonuniform_report as NR

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
    """
    Alternative module: Use the strategy in Bis-SNP to trim 5' bisulfite conversion failures
    """
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

    loc_dict = {}
    if filter_dup:
        ## if filter_up is TRUE, the duplicate reads will be assessed and shown in Dup_dis.pdf
        info("The filter_dup has been set True.")
        info("Assess the duplicate reads...")
        for sam in sam_inf:
            #check the input mapping files
            sam_format, read_inf = check.check_mapping_file(sam, s_path)
            if single_on:
                for read in read_inf:
                    loc_dict = LI.Loc_single(read, loc_dict, bsm)
            else:
                for read in read_inf:
                    loc_dict = LI.Loc_paired(read, loc_dict, bsm)
        max_cov = DR.duplicate_report(loc_dict, gsize, p_poisson, name)
        info('Get the duplicate reads distribution!')

    #get reference information
    ref = GR.get_ref(ref_file)
    trim_position = []

    filter_duplicate_reads = 0
    filter_nonuniform_trim_bp = 0
    filter_nonuniform_trim_bp_CG = 0
    filter_remove_overlap_bp = 0
    filter_not_mapping_reads = 0
    all_reads = 0
    not_mapping_reads = 0
    all_mapping_bp = 0

    ##filter the 5' bisulfite failure
    for sam in sam_inf:
        out_sam = sam[:-4] + '_' + name + '_filter.sam'
        out = open(out_sam, 'w')
        #check the input mapping files
        record_mate = {}
        sam_format, read_inf = check.check_mapping_file_header(sam, s_path)

        for read in read_inf:
            #for sam header
            if read.startswith('@'):
                out.write(read)
                continue
            else:
                all_reads += 1  ##record the read number (2013-06-20)

                #Get the read information for trimming
                #If the read isn't unique mapping, we will get a empty list ([]).
                #In: single unique mapping read  Out: [flag,strand,chr,pos,CIGAR,seq,score]
                #In: paired unique mapping read  Out: [flag,strand,chr,pos1,CIGAR,pos2,insert,seq,score]
            read_info = RI(read, bsm)
            read_info = read_info.extract_information()

            if len(read_info) == 0:
                not_mapping_reads += 1
                if not_mapping:         #keep the not_unique mapping reads (or not paired mapping)
                    out.write(read)
                else:
                    filter_not_mapping_reads += 1  ##record the not mapping read number (2013-06-20)
                continue

            if len(loc_dict) > 0: #the --filter_dup has been set True, have to remove duplicate reads
                duplicate, loc_dict = DF(read_info, loc_dict, max_cov, single_on)
            else:
                duplicate = False

            if single_on:
                all_mapping_bp += len(read_info[5])   ##record the mapping read basepair (2013-06-20)
            else:
                all_mapping_bp += len(read_info[7])   ##record the mapping read basepair (2013-06-20)

            record_mate, trim_position, filter_nonuniform_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp = NF.nonuniform_filter(read,
                                                                                                                out,
                                                                                                                read_info,
                                                                                                                ref,
                                                                                                                remove_overlap,
                                                                                                                duplicate,
                                                                                                                single_on,
                                                                                                                record_mate,
                                                                                                                trim_position,
                                                                                                                filter_nonuniform_trim_bp_CG,
                                                                                                                filter_duplicate_reads,
                                                                                                                filter_remove_overlap_bp)
        out.close()
        del record_mate
    NR.nonuniform_generator(trim_position, name)

    for i in range(len(trim_position)):
        filter_nonuniform_trim_bp += i * trim_position[i]

    ##produce the filter report
    info('Produce the report file...')
    report_out = open(name + "_BSeQC_nonuniform_filter_report.txt", 'w')
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
        report_out.write('Filter Duplicate reads: %d(%.2f%s of unique mapping reads)\n' % (
            filter_duplicate_reads, float(filter_duplicate_reads) / (all_reads - not_mapping_reads) * 100, "%"))
        report_out.write("Filter 5' nonconversion basepairs: %d(%.2f%s of unique mapping basepairs)\n" % (
            filter_nonuniform_trim_bp, float(filter_nonuniform_trim_bp) / all_mapping_bp * 100, "%"))
        report_out.write("Filter 5' nonconversion CpG basepairs: %d(%.2f%s of unique mapping basepairs)\n" % (
            filter_nonuniform_trim_bp_CG, float(filter_nonuniform_trim_bp_CG) / all_mapping_bp * 100, "%"))

    else:
        report_out.write('Not unique paired mapping reads: %d(%.2f%s)\n' % (
            not_mapping_reads, float(not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('Unique paired mapping reads: %d(%.2f%s)\n' % (
            (all_reads - not_mapping_reads), float(all_reads - not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('Skip not paired unique mapping reads: %d(%.2f%s)\n' % (
            filter_not_mapping_reads, float(filter_not_mapping_reads) / all_reads * 100, "%"))
        report_out.write('In unique paired mapping reads:\n')
        report_out.write('All unique paired mapping basepairs: %d\n' % all_mapping_bp)
        report_out.write('Filter Duplicate reads: %d(%.2f%s of unique paired mapping reads)\n' % (
            filter_duplicate_reads, float(filter_duplicate_reads) / (all_reads - not_mapping_reads * 100), "%"))
        report_out.write("Filter 5' nonconversion basepairs: %d(%.2f%s of unique mapping basepairs)\n" % (
            filter_nonuniform_trim_bp, float(filter_nonuniform_trim_bp) / all_mapping_bp * 100, "%"))
        report_out.write("Filter 5' nonconversion CpG basepairs: %d(%.2f%s of unique mapping basepairs)\n" % (
            filter_nonuniform_trim_bp_CG, float(filter_nonuniform_trim_bp_CG) / all_mapping_bp * 100, "%"))
        report_out.write('Filter overlapped basepairs: %d(%.2f%s of unique paired mapping basepairs)\n' % (
            filter_remove_overlap_bp, float(filter_remove_overlap_bp) / all_mapping_bp * 100, "%"))
    report_out.close()
    info('Get the report file!')