#!/usr/bin/python

# ------------------------------------
#python package
# ------------------------------------
import numpy as np

# ------------------------------------
#own python modules
# ------------------------------------
from BSeQC.read.read_info import MappingReader as RI


class MethReader:
    '''
    For a single-end or paired-end read, according to the reference genome,
    if the C in one position is CpG, we record the methylation state.
    the process is stratified by strand and read length.
    '''

    def __init__(self, read, strand_p, ref, bsmb, original_length, sequence_context):
        '''
        Initialize
        '''
        self.read = read
        self.strand_p = strand_p
        self.ref = ref
        self.bsmb = bsmb
        self.original_length = original_length
        self.sequence_context = sequence_context
        self.BS_conversion = {'+': ('C', 'T'), '-': ('G', 'A')}
        self.reverse_strand = ['-+']

    def get_read_info(self):
        '''
        Using the read class
        If the read is unique mapping or unique and  paired mapping, we will get a information list.
        If not, we will get a empty list ([]).
        In: single unique mapping read  Out: [flag,strand,chr,pos,CIGAR,seq,score]
        In: paired unique mapping read  Out: [flag,strand,chr,pos1,CIGAR,pos2,insert,seq,score]
        '''
        read_s = RI(self.read, self.bsmb)
        read_info = read_s.extract_information()
        return read_info


    def get_length(self, readlen, strand):
        '''
        If the user give the original length, we will use it, and ignore the read length.
        In some situation, The read length maybe not the sequence length, because of the the adapter cutting or other reasons
        '''
        if self.original_length:
            length = int(self.original_length)
        else:
            length = readlen
        return length

    def stratify_length(self, strand, length):
        '''
        stratify the read length to record the methylation state of CpG for each strand
        '''

        if not self.strand_p[strand].has_key(length):
            self.strand_p[strand][length] = np.zeros((2, int(length)))

        return self.strand_p

    def scan_context(self, seq, refseq, strand, length, match, convert, index):
        '''
        scan both the CG and the nonCG  (2013-06-04)
        '''
        if self.sequence_context == 'CG':
            while index >= 0:
                ref_index = index + 1
                if strand[1] == '+':                             # '++' or '-+'
                    context = refseq[ref_index:ref_index + 2]
                else:
                    context = refseq[ref_index - 1:ref_index + 1]     # '+-' or '+-'
                if (context == 'CG' and strand not in self.reverse_strand) or (
                            context == 'GC' and strand in self.reverse_strand):
                    if seq[index] == convert:
                        self.strand_p[strand][length][1, index] += 1
                    elif seq[index] == match:
                        self.strand_p[strand][length][1, index] += 1
                        self.strand_p[strand][length][0, index] += 1
                index = refseq[1:-1].find(match, index + 1)
        else:
            while index >= 0:
                ref_index = index + 1
                if strand[1] == '+':                                  # '++' or '-+'
                    context = refseq[ref_index:ref_index + 2]
                else:
                    context = refseq[ref_index - 1:ref_index + 1]     # '+-' or '+-'
                if (context != 'CG' and strand not in self.reverse_strand) or (
                            context != 'GC' and strand in self.reverse_strand):
                    if seq[index] == convert:
                        self.strand_p[strand][length][1, index] += 1
                    elif seq[index] == match:
                        self.strand_p[strand][length][1, index] += 1
                        self.strand_p[strand][length][0, index] += 1
                index = refseq[1:-1].find(match, index + 1)
        return self.strand_p

    def record_meth(self):
        '''
        get the read information to record the methylation state
        '''
        read_info = self.get_read_info()
        if len(read_info) == 0:
            return self.strand_p
        strand, chr, pos, seq = read_info[1], read_info[2], int(read_info[3]) - 1, read_info[-2]

        if len(chr.split('_')) > 1:   #filter some 'chr_' chromosome
            return self.strand_p

        readlen = len(seq)
        refseq = self.ref[chr][pos - 1:pos + readlen + 1]
        length = self.get_length(readlen, strand)
        self.strand_p = self.stratify_length(strand, length)

        # if the strand is '-+' (single-end) or '+-', reverse the seq and refseq
        if strand in self.reverse_strand:
            seq = seq[::-1]
            refseq = refseq[::-1]

        # record the methylation state for each unique (paired) mapping read

        match, convert = self.BS_conversion[strand[0]]
        index = refseq[1:-1].find(match) # for recognizing the CG and the non-CG site
        self.strand_p = self.scan_context(seq, refseq, strand, length, match, convert, index)
        return self.strand_p


class SingleMethReader(MethReader):
    '''
    specialise in single-end read
    '''

    def __init__(self, read, strand_p, ref, bsmb, original_length, sequence_context):
        MethReader.__init__(self, read, strand_p, ref, bsmb, original_length, sequence_context)


class PairMethReader(MethReader):
    '''
    specialise in paired-end read
    '''

    def __init__(self, read, strand_p, ref, bsmb, original_length, sequence_context):
        MethReader.__init__(self, read, strand_p, ref, bsmb, original_length, sequence_context)
        self.reverse_strand = ['-+', '+-']

    def get_length(self, readlen, strand):
        '''
        If the user give the original length, we will use it, and ignore the read length.
        The read length maybe not the sequence length, because of the the adapter cutting or other reasons
        '''
        if strand == '++' or strand == '-+':
            if self.original_length:
                length = int(self.original_length.split('_')[0])
            else:
                length = readlen
        else:
            if self.original_length:
                length = int(self.original_length.split('_')[1])
            else:
                length = readlen
        return length

    def record_meth(self):
        '''
        For paired-end read, if pos1 == pos2, the read will be ignored.
        '''
        read_info = self.get_read_info()
        if len(read_info) == 0:
            return self.strand_p
        #pos1 == pos2?
        elif int(read_info[3]) == int(read_info[5]):
            return self.strand_p
        else:
            MethReader.record_meth(self)
            return self.strand_p

