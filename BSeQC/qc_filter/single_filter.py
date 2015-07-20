#!/usr/bin/python


'''
trim the mbias in the single-end read
'''

def read_trim_single(read, strand_t, out, read_info, original_length, duplicate, filter_mbias_trim_bp,
                     filter_duplicate_reads):
    '''
    specialise in trimming single-end read
    '''
    if duplicate:#remove duplicate reads
        filter_duplicate_reads += 1   ##record the duplicate read (2013-06-20)
        return filter_mbias_trim_bp, filter_duplicate_reads

    strand, pos, CIGAR, seq, score = read_info[1], int(read_info[3]), read_info[4], read_info[5], read_info[6]
    if original_length:
        readlen = original_length
    else:
        readlen = len(seq)

    pos, seq, score, CIGAR = get_trmmed_positions(strand_t, strand, pos, CIGAR, seq, score, readlen)
    filter_mbias_trim_bp += len(read_info[5]) - len(seq)  ##record the mbias trim basepair (2013-06-20)

    #output the new read
    original_read_info = read.rstrip().split('\t')
    out.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s'
              % (original_read_info[0], original_read_info[1], original_read_info[2], pos, original_read_info[4],
                 CIGAR, original_read_info[6], original_read_info[7], original_read_info[8], seq, score))
    for i in range(11, len(original_read_info)):
        out.write('\t%s' % original_read_info[i])
    out.write('\n')
    return filter_mbias_trim_bp, filter_duplicate_reads


def get_trmmed_positions(strand_t, strand, pos, CIGAR, seq, score, readlen):
    #get the trimming positions
    #note: the length of some mapping reads may be shorter than the original sequence length
    if len(seq) >= (strand_t[strand][readlen][1] + 1):
        exclude_3 = strand_t[strand][readlen][1]
    else:
        exclude_3 = len(seq) - 1
    exclude_5 = strand_t[strand][readlen][0]

    if strand == '-+':
        pos = pos + (len(seq) - exclude_3 - 1)
        seq = seq[::-1][exclude_5:exclude_3 + 1][::-1]
        score = score[::-1][exclude_5:exclude_3 + 1][::-1]
        CIGAR = str(len(seq)) + CIGAR[-1]
    else:
        pos = pos + exclude_5
        seq = seq[exclude_5:exclude_3 + 1]
        score = score[exclude_5:exclude_3 + 1]
        CIGAR = str(len(seq)) + CIGAR[-1]

    return pos, seq, score, CIGAR
