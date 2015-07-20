#!/usr/bin/python



def read_trim_paired(read, ref, strand_t, out, read_info, original_length, auto, remove_overlap, duplicate,
                     filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp):
    '''
    trim the Mbias or only keep one copy of the overlapping segment in the paired-end read
    '''

    if duplicate:
        filter_duplicate_reads += 1   ##record the duplicate read (2013-06-20)
        return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

    strand, chr, pos1, CIGAR, pos2, insert, seq, score =\
    read_info[1], read_info[2], int(read_info[3]), read_info[4], int(read_info[5]), abs(int(read_info[6])), read_info[7], read_info[8]

    #only remove overlap, not trimming Mbias

    if not auto and remove_overlap:
        if strand == '++' or strand == '--':
            CIGAR, seq, score, filter_remove_overlap_bp = remove_overlap_seq(pos1, pos2, insert, seq, score, CIGAR,
                filter_remove_overlap_bp)
            if len(seq) == 0:
            #In this situation, the two read mates are the same, we should count one time if remove_overlap is True
                return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp
            else:
                writer(pos1, CIGAR, pos2, insert, seq, score, read, out)
                return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp
        else:
            out.write(read)
            return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp


    #get the readlen and the corresponding mate length
    readlen, p_len = get_p_len(original_length, strand, pos1, pos2, seq, insert)

    if p_len <= 0:   #2014.08.12
        out.write(read)
        return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

    if not strand_t[strand].has_key(readlen) and (pos1 == pos2):
        if remove_overlap and (strand == '++' or strand == '--'):
            filter_remove_overlap_bp += len(read_info[7])
            return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp
        else:
            out.write(read)
            return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

    if (pos2 < pos1 and (strand == '++' or strand == '--')) or (pos1 < pos2 and (strand == '+-' or strand == '-+')):
        out.write(read)
        return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp

    #get the trimming positions
    #using the other paired read to decide the new pos2 and insert size
    pos1_new, pos2_new, CIGAR, insert, seq, score, filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_remove_overlap_bp = get_trimmed_position(
        strand_t, ref, strand, chr, readlen, p_len,
        pos1, pos2, seq, score, CIGAR, remove_overlap, filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_remove_overlap_bp)
    if len(seq) == 0:
        return filter_mbias_trim_bp,filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp
    else:
        writer(pos1_new, CIGAR, pos2_new, insert, seq, score, read, out)
    return filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_duplicate_reads, filter_remove_overlap_bp


def remove_overlap_seq(pos1, pos2, insert, seq, score, CIGAR, filter_remover_overlap_bp):
    ## If two reads mate are overlapped with each other, only keep one copy of the the overlapping segment
    if insert < len(seq) + (insert - (pos2 - pos1)):
        CIGAR = str(pos2 - pos1) + CIGAR[-1]
        filter_remover_overlap_bp += len(seq) - len(seq[:(pos2 - pos1)]) ##record the overlap basepair (2013-06-20)
        seq = seq[:(pos2 - pos1)]
        score = score[:(pos2 - pos1)]
    return CIGAR, seq, score, filter_remover_overlap_bp


def writer(pos1, CIGAR, pos2, insert, seq, score, read, out):
    '''
    output the trimmed or overlap-removed read
    '''
    original_read_info = read.rstrip().split('\t')
    out.write('%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%s'
              % (original_read_info[0], original_read_info[1], original_read_info[2], pos1, original_read_info[4],
                 CIGAR, original_read_info[6], pos2, insert, seq, score))
    for i in range(11, len(original_read_info)):
        out.write('\t%s' % original_read_info[i])
    out.write('\n')
    return


def get_p_len(original_length, strand, pos1, pos2, seq, insert):
    '''
    get the readlen and the corresponding mate length
    '''
    if original_length:
        if strand == '++' or strand == '-+':
            readlen = original_length[0]
            p_len = original_length[1]
        else:
            readlen = original_length[1]
            p_len = original_length[0]
    else:
        readlen = len(seq)
        p_len = insert - abs(pos2 - pos1)
    return readlen, p_len


def get_trimmed_position(strand_t, ref, strand, chr, readlen, p_len, pos1, pos2, seq, score, CIGAR, remove_overlap,
                         filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_remove_overlap_bp):
    #get the trimming positions
    #using the other paired read to decide the new pos2 and insert size

    p_dict = {'++': '+-', '+-': '++', '-+': '--', '--': '-+'}
    p_strand = p_dict[strand]
    p_exclude_5 = strand_t[p_strand][p_len][0]
    p_exclude_3 = strand_t[p_strand][p_len][1]
    exclude_5 = strand_t[strand][readlen][0]
    exclude_3 = strand_t[strand][readlen][1]

    if strand == '++' or strand == '--':
        pos1_new = pos1 + exclude_5
        pos2_new = pos2 + (p_len - 1 - p_exclude_3)  #the correspond pair read
        pos1_3_new = pos1 + exclude_3
        pos2_3_new = pos2 + (p_len - 1 - p_exclude_5)
        #after trimming if 5' > 3'
        if pos1_new > pos2_new:
            pos2_new = pos1_new
        if pos1_3_new > pos2_3_new:
            pos1_3_new = pos2_3_new
        start = pos1_new - pos1
        end = pos1_3_new - pos1 + 1
        filter_mbias_trim_bp += len(seq) - len(seq[start:end])  ##record the mbias basepair (2013-06-20)
        ref_seq = ref[chr][pos1 - 2: pos1 - 1 + readlen]
        ref_seq_trim = ref_seq[:start + 2]
        filter_mbias_trim_bp_CG = get_trim_CG(ref_seq_trim, strand, filter_mbias_trim_bp_CG)
        seq = seq[start:end]
        score = score[start:end]
        CIGAR = str(len(seq)) + CIGAR[-1]
        #add the mate read length
        insert = pos2_new - pos1_new + (pos2_3_new - pos2_new + 1)

        ## If two reads mate are overlapped with each other, only keep one copy of the the overlapping segment
        if remove_overlap:
            CIGAR, seq, score, filter_remove_overlap_bp = remove_overlap_seq(pos1_new, pos2_new, insert, seq, score,
                CIGAR, filter_remove_overlap_bp)
    else:
        pos1_new = pos1 + (readlen - 1 - exclude_3)
        pos2_new = pos2 + p_exclude_5
        pos1_3_new = pos1 + (readlen - 1 - exclude_5)
        #pos2_3_new = pos2 + p_exclude_3
        #after trimming if 5' > 3'
        if pos1_new < pos2_new:
            pos1_new = pos2_new
            #if pos1_3_new < pos2_3_new:
            #pos2_3_new = pos1_3_new
        start = pos1_new - pos1
        end = pos1_3_new - pos1 + 1
        filter_mbias_trim_bp += len(seq) - len(seq[start:end])  ##record the mbias basepair (2013-06-20)
        ref_seq = ref[chr][pos1 - 2: pos1 + readlen]
        ref_seq_trim = ref_seq[end:]
        filter_mbias_trim_bp_CG = get_trim_CG(ref_seq_trim, strand, filter_mbias_trim_bp_CG)
        seq = seq[start:end]
        score = score[start:end]
        CIGAR = str(len(seq)) + CIGAR[-1]
        #add self read length
        insert = -(pos1_new - pos2_new + (pos1_3_new - pos1_new + 1))

    return pos1_new, pos2_new, CIGAR, insert, seq, score, filter_mbias_trim_bp, filter_mbias_trim_bp_CG, filter_remove_overlap_bp

def get_trim_CG(ref_seq_trim, strand, filter_mbias_trim_bp_CG):
    BS_conversion = {'+': ('C', 'T'), '-': ('G', 'A')}
    reverse_strand = ['-+', '+-']
    match, convert = BS_conversion[strand[0]]
    index = ref_seq_trim[1:-1].find(match)
    while index >= 0:
        ref_index = index + 1
        if strand[0] == '+':                             # '++' or '+-'
            context = ref_seq_trim[ref_index:ref_index + 2]
        else:
            context = ref_seq_trim[ref_index - 1:ref_index + 1]

        if context == 'CG':
            filter_mbias_trim_bp_CG += 1
        index = ref_seq_trim[1:-1].find(match, index + 1)
    return filter_mbias_trim_bp_CG
