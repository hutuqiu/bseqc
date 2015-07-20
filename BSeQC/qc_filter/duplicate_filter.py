#!/usr/bin/python



def duplicate_filter(read_info,loc_dict,max_cov,single_on):
    '''
    filter the duplicate reads, let the clonal reads not more than max_cov
    return:
    duplicate: True or False
    loc_dict: If the read is duplicate, the correlated site in loc_dict minus the filter_shift
              until the coverage is not more than max_cov
              (1 for single-end, 0.5 for paired-end).
    '''
    if single_on:
        filter_shift = 1
        strand, chrom, pos, seq= read_info[1], read_info[2], read_info[3], read_info[5]
        if strand == '-+':
            pos = int(pos) + len(seq) -1
        site = strand + '_' + str(pos)
    else:
        filter_shift = 0.5
        strand, chrom, pos, pos2, insert = read_info[1], read_info[2], read_info[3], read_info[5], read_info[6]
        if strand == '++' or strand == '--':
            site = pos + '_' + str(int(pos) + int(insert) - 1)
        else:  #for '+-' or '-+' strand, the insert size is negative
            site = pos2 + '_' + str(int(pos2) + abs(int(insert)) - 1)

    ##for some strange mapping
    if (strand == '-+' or strand == '+-') and not loc_dict[chrom].has_key(site):
        duplicate = False
        return duplicate, loc_dict


    if loc_dict[chrom][site] > max_cov:
        duplicate = True
        loc_dict[chrom][site] -= filter_shift
    else:
        duplicate = False
    return duplicate, loc_dict