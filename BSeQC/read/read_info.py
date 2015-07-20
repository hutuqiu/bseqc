#!/usr/bin/python


class MappingReader:
    '''Get the information from mapping read

    '''

    def __init__(self, read, bsmb):
        '''
        Initialize the read information
        '''
        self.read_info = read.rstrip().split('\t')
        self.bsmb = bsmb

    def flag_conversion(self):
        '''
        convert the decimal flag to a string for reading
        example:
        83 -> pPr1
        99 -> pPR1
        '''
        flag_d = self.read_info[1]                                          # the decimal flag in sam file
        flag_b = bin(int(flag_d))[2:][::-1]                                 # convert a decimal flag to a binary flag
        flag_c = ['p', 'P', 'u', 'U', 'r', 'R', '1', '2', 's', 'f', 'd']              # for a binary flag conversion
        flag = [flag_c[i] for i in range(len(flag_b)) if flag_b[i] == '1']    # convert a binary flag to a string
        return flag

    def extract_strand(self):
        '''
        using the flag to extract the mapping strand information
        '''

        if self.bsmb:
            if self.read_info[-1].split(':')[-1] == 'W':
                strand = '++'
            elif self.read_info[-1].split(':')[-1] == 'WR':
                strand = '+-'
            elif self.read_info[-1].split(':')[-1] == 'C':
                strand = '-+'
            elif self.read_info[-1].split(':')[-1] == 'CR':
                strand = '--'

        else:
            flag = self.flag_conversion()
            if '1' in flag:                     #for the paired end
                #if 'R' in flag:
                #   strand = '++'
                if 'r' in flag:                  #support Bismark mapping result (2013-06-02)
                    strand = '-+'
                else:
                    strand = '++'
            elif '2' in flag:                    #for the paired end
                if 'R' in flag:
                    strand = '--'
                else:
                    strand = '+-'
            else:                                 #for single end
                if 'r' in flag:
                    strand = '-+'
                else:
                    strand = '++'
        return strand

    def extract_information(self):
        '''
        If the read is unique mapping or unique and  paired mapping, we will get a information list.
        If not, we will get a empty list ([]).
        In: single unique mapping read  Out: [flag,strand,chr,pos,CIGAR,seq,score]
        In: paired unique mapping read  Out: [flag,strand,chr,pos1,CIGAR,pos2,insert,seq,score]
        '''
        flag = self.flag_conversion()
        read_info = []
        if 'u' in flag or 's' in flag or 'd' in flag or 'f' in flag:       #add or 'd' 'f' in flag on Aug 19
            return read_info
        if 'p' in flag and 'P' not in flag:
            return read_info
        chr = self.read_info[2]
        strand = self.extract_strand()
        CIGAR = self.read_info[5]
        seq = self.read_info[9]
        score = self.read_info[10]
        if 'p' in flag:
            pos1 = self.read_info[3]
            pos2 = self.read_info[7]
            insert = self.read_info[8]
            read_info = [flag, strand, chr, pos1, CIGAR, pos2, insert, seq, score]
        else:
            pos = self.read_info[3]
            read_info = [flag, strand, chr, pos, CIGAR, seq, score]
        return read_info
