#!/usr/bin/python

# ------------------------------------
#python package
# ------------------------------------
import logging
import sys

# ------------------------------------
#the read package
# ------------------------------------
from BSeQC.read.read_info import MappingReader as RI
from BSeQC.qc_filter.single_filter import read_trim_single as SF
from BSeQC.qc_filter.paired_filter import read_trim_paired as PF
from BSeQC.qc_filter.duplicate_filter import duplicate_filter as DF
from BSeQC.read import check_file as check
from BSeQC.read import get_reference as GR


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


def filter_sam(sam_inf, ref_file, bsmb, strand_t, read_l, single_on, name, s_path, auto, remove_overlap, loc_dict, max_cov,
               not_mapping):
    '''
    Trim the mapping files with the biased positions of every length in every strand,
    which are saved in the variance: strand_t.
    '''
    filter_duplicate_reads = 0
    filter_mbias_trim_bp = 0
    filter_mbias_trim_bp_CG = 0
    filter_remove_overlap_bp = 0
    filter_not_mapping_reads = 0
    all_reads = 0
    not_mapping_reads = 0
    all_mapping_bp = 0
    ref = GR.get_ref(ref_file)
    for sam in sam_inf:
        out_sam = sam[:-4] + '_' + name + '_filter.sam'
        out = open(out_sam, 'w')

        #check the input mapping files
        sam_format, read_inf = check.check_mapping_file_header(sam, s_path)

        #scan every read to qc_filter
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
            read_info = RI(read, bsmb)
            read_info = read_info.extract_information()
            if len(read_info) == 0:
                not_mapping_reads += 1
                if not_mapping:         #keep the not_unique mapping reads
                    out.write(read)
                else:
                    filter_not_mapping_reads += 1  ##record the not mapping read number (2013-06-20)
                continue

            if len(loc_dict) > 0: #the --filter_dup has been set True, have to remove duplicate reads
                duplicate, loc_dict = DF(read_info, loc_dict, max_cov, single_on)
            else:
                duplicate = False

            if single_on:
                all_mapping_bp += len(read_info[5])     ##record the mapping read basepair (2013-06-20)
                if auto:
                    if read_l[0] != '':
                        original_length = int(read_l[sam_inf.index(sam)])
                    else:
                        original_length = ''
                    filter_mbias_trim_bp, filter_duplicate_reads = SF(read, strand_t, out, read_info, original_length,
                        duplicate, filter_mbias_trim_bp, filter_duplicate_reads)
                else:
                    if not duplicate and len(loc_dict) > 0:
                        out.write(read)                 #not trimming, only output not_duplicate reads
                    else:
                        filter_duplicate_reads += 1     ##record the duplicate read (2013-06-20)
            else:
                all_mapping_bp += len(read_info[7])     ##record the mapping read basepair (2013-06-20)
                if auto or remove_overlap:
                    if read_l[0] != '':
                        original_length = [int(i) for i in read_l[sam_inf.index(sam)].split('_')]
                    else:
                        original_length = ''
                    filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp = PF(read, ref, strand_t, out,
                        read_info, original_length, auto, remove_overlap, duplicate, filter_mbias_trim_bp, filter_mbias_trim_bp_CG,
                        filter_duplicate_reads, filter_remove_overlap_bp)
                else:
                    if not duplicate and len(loc_dict) > 0:
                        out.write(read)                  #not trimming, only output not_duplicate reads
                    else:
                        filter_duplicate_reads += 1     ##record the duplicate read (2013-06-20)
        out.close()

    ##produce the filter report
    info('Produce the report file...')
    report_out = open(name + "_BSeQC_mbias_filter_report.txt", 'w')
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
        #report_out.write('Filter Mbias CpG basepairs: %d(%.2f%s of unique mapping basepairs)\n' % (
        #filter_mbias_trim_bp_CG, float(filter_mbias_trim_bp_CG) / all_mapping_bp * 100, "%"))
        report_out.write('Filter Mbias basepairs: %d(%.2f%s of unique mapping basepairs)\n' % (
        filter_mbias_trim_bp, float(filter_mbias_trim_bp) / all_mapping_bp * 100, "%"))

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
        filter_duplicate_reads, float(filter_duplicate_reads) / (all_reads - not_mapping_reads) * 100, "%"))
        report_out.write('Filter Mbias basepairs: %d(%.2f%s of unique paired mapping basepairs)\n' % (
        filter_mbias_trim_bp, float(filter_mbias_trim_bp) / all_mapping_bp * 100, "%"))
        report_out.write("Filter 5' Mbias CpG basepairs: %d(%.2f%s of unique mapping basepairs)\n" % (
        filter_mbias_trim_bp_CG, float(filter_mbias_trim_bp_CG) / all_mapping_bp * 100, "%"))
        report_out.write('Filter overlapped basepairs: %d(%.2f%s of unique paired mapping basepairs)\n' % (
        filter_remove_overlap_bp, float(filter_remove_overlap_bp) / all_mapping_bp * 100, "%"))
    report_out.close()
    info('Get the report file!')




	
