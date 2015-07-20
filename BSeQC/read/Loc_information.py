#!/usr/bin/python

'''
Record the location information for removing duplicate reads
'''


# ------------------------------------
#own python modules
# ------------------------------------
from BSeQC.read.read_info import MappingReader as RI

def Loc_single(read,loc_dict,bsmb):
    '''
    Record the location information for single-end sequencing
    Return:
    loc_s = { chrid: {loc1: coverage; ...}; chrid:..}
    loc1: 5' pos + strand
    '''

    # Using the read class  In: single unique mapping read  Out: [flag,strand,chrom,pos,CIGAR,seq,score]
    read_s = RI(read, bsmb)
    read_info = read_s.extract_information()
    if len(read_info) == 0:
        return loc_dict

    strand, chrom, pos, seq = read_info[1], read_info[2], read_info[3], read_info[5]

    # The minus strand should be shifted len(seq) bp
    if strand == '-+':
        pos = int(pos) +  len(seq) -1

    site = strand + '_' + str(pos)
    if loc_dict.has_key(chrom):
        if loc_dict[chrom].has_key(site):
            loc_dict[chrom][site] += 1
        else:
            loc_dict[chrom][site] = 1
    else:
        loc_dict[chrom] = {}
        loc_dict[chrom][site] = 1
    return loc_dict

def Loc_paired(read,loc_dict,bsmb):
    '''
    Record the location information for paired sequencing
    Return:
    loc_s = { chrid: {loc1: coverage; ...}; chrid:..}
    loc1: plus strand (++ or --) 5' pos +  minus strand (+- or -+) 5' pos
    '''

    # Using the read class In: paired unique mapping read  Out: [flag,strand,chrom,pos1,CIGAR,pos2,insert,seq,score]
    read_s = RI(read,bsmb)
    read_info = read_s.extract_information()
    if len(read_info) == 0:
        return loc_dict
    strand, chrom, pos, pos2, insert = read_info[1], read_info[2], read_info[3], read_info[5], read_info[6]


    # only use one mate with plus insert to record
    if strand == '+-' or strand == '-+':
        return loc_dict

    site = pos + '_' + str(int(pos) + int(insert) - 1)
    if loc_dict.has_key(chrom):
        if loc_dict[chrom].has_key(site):
            loc_dict[chrom][site] +=1
        else:
            loc_dict[chrom][site] = 1
    else:
        loc_dict[chrom] = {}
        loc_dict[chrom][site] = 1
    return loc_dict
