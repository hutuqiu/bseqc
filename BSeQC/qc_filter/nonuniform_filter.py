# ------------------------------------
# python modules
# ------------------------------------
import os
import sys
import logging

# ------------------------------------
#own python modules
# ------------------------------------
from BSeQC.qc_filter import paired_filter as PF



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


def nonuniform_filter(read, out, read_info, ref, remove_overlap, duplicate, single_on, record_mate, trim_position, filter_nonuniform_trim_bp_CG,
                      filter_duplicate_reads, filter_remove_overlap_bp):
    '''
    Use the strategy in Bis-SNP to trim 5' bisulfite conversion failures
    For each read, we walk along the read from 5' to 3',
    and we remove any Cs on the C-strand until we reach the first reference C which is converted to a T.
    '''
    #if the read is duplicate read and the --filter_dup has been set True, the read will be removed
    if duplicate:
        filter_duplicate_reads += 1   ##record the duplicate read (2013-06-20)
        return record_mate, trim_position, filter_nonuniform_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

    if single_on:
        strand, chr, pos1, CIGAR, seq, score = read_info[1], read_info[2], int(read_info[3]) - 1, read_info[4], \
                                               read_info[5], read_info[6]
    else:
        strand, chr, pos1, CIGAR, pos2, insert, seq, score = \
            read_info[1], read_info[2], int(read_info[3]) - 1, read_info[4], int(read_info[5]) - 1, abs(
                int(read_info[6])), read_info[7], read_info[8]
        if pos2 < pos1 and (strand == '++' or strand == '--') or (pos1 < pos2 and (strand == '+-' or strand == '-+')):
            out.write(read)
            return record_mate, trim_position, filter_nonuniform_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

    BS_conversion = {'+': ('C', 'T'), '-': ('G', 'A')}
    reverse_strand = ['-+', '+-']
    original_read_info = read.rstrip().split('\t')

    #filter 5' conversion failure read by read
    readlen = len(seq)
    refseq = ref[chr][pos1:pos1 + readlen]
    refseq_for_trim = ref[chr][pos1 - 1:pos1 + readlen + 1]


    # if the strand is '-+' (single-end) or '+-', reverse the seq and refseq
    if strand in reverse_strand:
        seq = seq[::-1]
        refseq = refseq[::-1]
        score = score[::-1]

    match, convert = BS_conversion[strand[0]]
    index = refseq.find(match) #record the position

    while index >= 0:
        if seq[index] == convert:
            break
        else:
            index = refseq.find(match, index + 1)

    if index < 0: #no finding C or G
        index = 0

    if index < len(trim_position):
        trim_position[index] += 1
    else:
        new_trim_position = [0] * ((index + 1) - len(trim_position))
        trim_position = trim_position + new_trim_position
        trim_position[index] += 1


    #trimming
    if strand in reverse_strand:
        seq = seq[index:][::-1]
        score = score[index:][::-1]
        pos1 += 1
        if not single_on:
            pos2 += 1
    else:
        seq = seq[index:]
        score = score[index:]
        pos1 = pos1 + index + 1
        if not single_on:
            pos2 += 1
    CIGAR = str(len(seq)) + CIGAR[-1]


    if strand in reverse_strand:
        refseq_trim = refseq_for_trim[len(seq):]
    else:
        refseq_trim = refseq_for_trim[:index + 2]

    filter_nonuniform_trim_bp_CG = PF.get_trim_CG(refseq_trim, strand, filter_nonuniform_trim_bp_CG)




    if single_on:
        out.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s'
                  % (original_read_info[0], original_read_info[1], original_read_info[2], pos1, original_read_info[4],
                     CIGAR, original_read_info[6], original_read_info[7], original_read_info[8], seq, score))
        for i in range(11, len(original_read_info)):
            out.write('\t%s' % original_read_info[i])
        out.write('\n')
    else:
        #get the original information

        strand, chr, pos1_ori, pos2_ori, insert = \
            read_info[1], read_info[2], int(read_info[3]) - 1, int(read_info[5]) - 1, abs(int(read_info[6]))

        if strand in reverse_strand:
            site_information = [chr, str(pos2_ori), str(pos1_ori), str(insert)]  #change pos2 and pos1 in '+-' or '-+'
        else:
            site_information = [chr, str(pos1_ori), str(pos2_ori), str(insert)]
        read_mate_site = '_'.join(site_information)
        p_dict = {'++': '+-', '+-': '++', '-+': '--', '--': '-+'}

        if not record_mate.has_key(read_mate_site):
            record_mate[read_mate_site] = [[read, strand, CIGAR, seq, score, pos1, pos2]]
        else:
            record_mate[read_mate_site].append([read, strand, CIGAR, seq, score, pos1, pos2])
            if len(record_mate[read_mate_site]) == 2 and p_dict[record_mate[read_mate_site][0][1]] == \
                    record_mate[read_mate_site][1][1]:
                #get the new position
                for mate in record_mate[read_mate_site]:
                    if mate[1] not in reverse_strand:
                        pos1_new = mate[5]
                        pos2_new = mate[6]
                    #get new insert
                for mate in record_mate[read_mate_site]:
                    if mate[1] in reverse_strand:
                        insert = pos2_new - pos1_new + len(mate[3])

                for mate in record_mate[read_mate_site]:
                    read, strand, CIGAR, seq, score = mate[0], mate[1], mate[2], mate[3], mate[4]
                    if remove_overlap:
                        if strand not in reverse_strand:
                            CIGAR, seq, score, filter_remove_overlap_bp = PF.remove_overlap_seq(pos1_new, pos2_new,
                                                                                                insert, seq, score,
                                                                                                CIGAR,
                                                                                                filter_remove_overlap_bp)
                            if len(seq) == 0:
                            #In this situation, the two read mates are the same, we should count one time if remove_overlap is True
                                continue
                            else:
                                PF.writer(pos1_new, CIGAR, pos2_new, insert, seq, score, read, out)
                        else:
                            PF.writer(pos2_new, CIGAR, pos1_new, -insert, seq, score, read, out)
                    else:
                        if strand not in reverse_strand:
                            PF.writer(pos1_new, CIGAR, pos2_new, insert, seq, score, read, out)
                        else:
                            PF.writer(pos2_new, CIGAR, pos1_new, -insert, seq, score, read, out)
                del record_mate[read_mate_site]

    return record_mate, trim_position, filter_nonuniform_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

