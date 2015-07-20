#!/usr/bin/python
# ------------------------------------
#python package
# ------------------------------------
import sys
import logging
import copy
import numpy as np
from scipy.stats import norm

# ------------------------------------
#own python modules
# ------------------------------------
from BSeQC.qc_assess import mbias_report as MR
from BSeQC.qc_assess import duplicate_report as DR
from BSeQC.read import C_information as CI
from BSeQC.read import Loc_information as LI
from BSeQC.read import get_reference as GR
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


class QC_Report:
    '''
    The main QC_Report: Mbias QC report
    '''

    def __int__(self, sam_inf, ref_file, bsmb, s_path, name, read_l, single_on, pvalue, drift, trim_file):
        '''
        Initialize the read information
        '''
        self.sam_inf = sam_inf
        self.ref_file = ref_file
        self.bsmb = bsmb
        self.s_path = s_path
        self.name = name
        self.read_l = read_l
        self.single_on = single_on
        self.pvalue = pvalue
        self.drift = drift
        self.trim_file = trim_file


    def decide_trim_bp(self, strand_p):
        '''
        determine how much bp to qc_filter
        using two cutoffs: pvalue < = 0.01; mean +- drift (default:2%)
        '''
        strand_t = {}
        for s in strand_p.keys():
            strand_t[s] = {}
            for l in strand_p[s].keys():
                #some read length only include a few reads, which are not enough to measure the Mbias
                if 0 in strand_p[s][l][1,]:
                    zero_index = np.where(strand_p[s][l][1,] == 0)
                    strand_p[s][l][1,][zero_index] = strand_p[s][l][1,][zero_index] + 1
                    #del strand_p[s][l]
                    #strand_t[s][l] = 'NA'
                    #continue
                strand_t[s][l] = [0, 0]
                Mcall = strand_p[s][l][0,] / strand_p[s][l][1,] * 100 # M/(M+U)*100%
                p_len = len(Mcall)
                p_standard = Mcall[round(p_len * 0.3):round(p_len * 0.7)]
                mean = p_standard.mean()
                std_v = p_standard.std()

                cutoff1 = self.pvalue
                cutoff2 = self.drift
                if std_v == 0:
                    strand_t[s][l] = [0, p_len - 1]
                    continue
                for pos_s in range(p_len):
                    if norm(loc=mean, scale=std_v).cdf(Mcall[pos_s]) <= 0.5:
                        CDF = norm(loc=mean, scale=std_v).cdf(Mcall[pos_s])
                    else:
                        CDF = 1 - norm(loc=mean, scale=std_v).cdf(Mcall[pos_s])
                    if CDF >= cutoff1 or abs(Mcall[pos_s] - mean) <= cutoff2:
                        strand_t[s][l][0] = pos_s
                        break
                end_index_list = range(p_len)
                end_index_list.reverse()
                for pos_e in end_index_list:
                    if norm(loc=mean, scale=std_v).cdf(Mcall[pos_e]) <= 0.5:
                        CDF = norm(loc=mean, scale=std_v).cdf(Mcall[pos_e])
                    else:
                        CDF = 1 - norm(loc=mean, scale=std_v).cdf(Mcall[pos_e])
                        #if abs(Mcall[pos_e] - mean) <= cutoff1 or abs(Mcall[pos_e] - mean)<=cutoff2:
                    if CDF >= cutoff1 or abs(Mcall[pos_e] - mean) <= cutoff2:
                        strand_t[s][l][1] = pos_e
                        break
        return strand_t

    def decide_final_trimming(self, strand_t_raw):
        '''
        decide final trimming from CG Mbias-plot and nonCG Mbias plot (2013-06-04)
        '''
        strand_t_final = {}
        for s in strand_t_raw[0].keys():
            strand_t_final[s] = {}
            for l in strand_t_raw[0][s].keys():
                strand_t_final[s][l] = [0, 0]
                if strand_t_raw[0][s][l][0] >= strand_t_raw[1][s][l][0]:
                    strand_t_final[s][l][0] = strand_t_raw[0][s][l][0]
                else:
                    strand_t_final[s][l][0] = strand_t_raw[1][s][l][0]

                if strand_t_raw[0][s][l][1] <= strand_t_raw[1][s][l][1]:
                    strand_t_final[s][l][1] = strand_t_raw[0][s][l][1]
                else:
                    strand_t_final[s][l][1] = strand_t_raw[1][s][l][1]
        return strand_t_final


    def parser_sambam(self, strand_p, ref):
        #Scan each SAM file to calculate the M fractions

        #add non-CpG Mbias plot (2013-06-04)
        strand_p_CG = copy.deepcopy(strand_p)
        strand_p_nonCG = copy.deepcopy(strand_p)
        for s in range(len(self.sam_inf)):
            sam_format, read_inf = check.check_mapping_file(self.sam_inf[s], self.s_path)
            #If the user give the original length, we will use it, and ignore the read length.
            #The read length maybe not the sequence length, because of the the adapter cutting or other reasons
            #you can set it by using the -l or --len option
            if self.read_l[0] != '':
                original_length = self.read_l[s]
            else:
                original_length = ''
            if self.single_on:
                for read in read_inf:
                    rm_single_CG = CI.SingleMethReader(read, strand_p_CG, ref, self.bsmb, original_length, 'CG')
                    rm_single_nonCG = CI.SingleMethReader(read, strand_p_nonCG, ref, self.bsmb,  original_length, 'nonCG')
                    strand_p_CG = rm_single_CG.record_meth()
                    strand_p_nonCG = rm_single_nonCG.record_meth()
            else:
                for read in read_inf:
                    rm_pair_CG = CI.PairMethReader(read, strand_p_CG, ref, self.bsmb, original_length, 'CG')
                    rm_pair_nonCG = CI.PairMethReader(read, strand_p_nonCG, ref, self.bsmb, original_length, 'nonCG')
                    strand_p_CG = rm_pair_CG.record_meth()
                    strand_p_nonCG = rm_pair_nonCG.record_meth()
        strand_p = [strand_p_CG, strand_p_nonCG]
        return strand_p

    def user_defined_trimming(self):
        strand_t = {}
        for line in open(self.trim_file):
            each = line.rstrip().split('\t')
            strand = each[0]
            length = int(each[1])
            trim_pos = [int(each[2]) - 1, int(each[3]) - 1]
            if strand_t.has_key(strand):
                strand_t[strand][length] = trim_pos
            else:
                strand_t[strand] = {}
                strand_t[strand][length] = trim_pos
        return strand_t

    def produce_final_trim_file(self, strand_t):
        trim_file = open(self.name + '_final_trim_file.txt', 'w')
        for s in strand_t.keys():
            for l in strand_t[s].keys():
                trim_file.write('%s\t%d\t%d\t%d\n' % (s, l, strand_t[s][l][0] + 1, strand_t[s][l][1] + 1))
        trim_file.close()
        return


    def generator(self):
        '''
        The main method in the class
        1. Draw Mbias plot and generate Mbias table fore per read length in per strand
        2. Decide the trimming positions based on the Mbias plot and generate the trimming file
        '''
        if len(self.trim_file) != 0:
            info("Used the trimming file from the user defined!! Ignore the step of automatically deciding trimming.")
            return self.user_defined_trimming()
        ref = GR.get_ref(self.ref_file)
        info("Calculate the M fraction for every position...")
        # check: are the input SAM files paired-end or single-end
        strand_p = {}
        if self.single_on:
            strand_p['++'] = {}
            strand_p['-+'] = {}
        else:
            strand_p['++'] = {}
            strand_p['-+'] = {}
            strand_p['+-'] = {}
            strand_p['--'] = {}
        strand_p = self.parser_sambam(strand_p, ref)
        #modify in 2013-06-04
        strand_t_raw = []
        name_context = [self.name + '_CG', self.name + '_nonCG']
        for i in range(len(strand_p)):
            strand_t_each = self.decide_trim_bp(strand_p[i])
            strand_t_raw.append(strand_t_each)
            MR.mbias_generator(strand_p[i], strand_t_each, name_context[i])
        strand_t = self.decide_final_trimming(strand_t_raw)
        self.produce_final_trim_file(strand_t)
        return strand_t


class QC_Report_Mias(QC_Report):
    '''
    Mbias QC report:
    1. Draw Mbias plot and generate Mbias table fore per read length in per strand
    2. Decide the trimming positions based on the Mbias plot and generate the trimming file
    '''

    def __init__(self, sam_inf, ref_file, bsmb, s_path, name, read_l, single_on, pvalue, drift, trim_file):
        QC_Report.__int__(self, sam_inf, ref_file, bsmb, s_path, name, read_l, single_on, pvalue, drift, trim_file)


class QC_Report_Mbias_Dup(QC_Report):
    '''
    The filter_dup has been set True.
    The QC report includes both Mibas and Dup:
    1. Draw Mbias plot and generate Mbias table fore per read length in per strand
    2. Decide the trimming positions based on the Mbias plot and generate the trimming file
    3. Show the duplicate reads distribution
    '''

    def __init__(self, sam_inf, ref_file, bsmb, s_path, name, read_l, single_on, pvalue, drift, trim_file, p_poisson, gsize):
        QC_Report.__int__(self, sam_inf, ref_file, bsmb, s_path, name, read_l, single_on, pvalue, drift, trim_file)
        self.p_poisson = p_poisson
        self.gsize = gsize

    def parser_sambam(self, strand_p, ref):
        #Scan each SAM file to calculate the M fractions and collection the location information to find duplicate reads
        loc_dict = {}
        strand_p_CG = copy.deepcopy(strand_p)
        strand_p_nonCG = copy.deepcopy(strand_p)
        for s in range(len(self.sam_inf)):
            sam_format, read_inf = check.check_mapping_file(self.sam_inf[s], self.s_path)
            #If the user give the original length, we will use it, and ignore the read length.
            #The read length maybe not the sequence length, because of the the adapter cutting or other reasons
            #you can set it by using the -l or --len option
            if self.read_l[0] != '':
                original_length = self.read_l[s]
            else:
                original_length = ''
                ## If the trim_file is used, the assessment of Mbias will be skipped
            if len(self.trim_file) != 0:
                info("Only record the location information of %s..." % self.sam_inf[s])
                if self.single_on:
                    for read in read_inf:
                        loc_dict = LI.Loc_single(read, loc_dict, self.bsmb)
                else:
                    for read in read_inf:
                        loc_dict = LI.Loc_paired(read, loc_dict, self.bsmb)
            else:
                #add non-CpG Mbias plot (2013-06-04)
                info("Calculate the M fraction for every position and record the location information in %s..." %
                     self.sam_inf[s])
                if self.single_on:
                    for read in read_inf:
                        rm_single_CG = CI.SingleMethReader(read, strand_p_CG, ref, self.bsmb, original_length, 'CG')
                        rm_single_nonCG = CI.SingleMethReader(read, strand_p_nonCG, ref, self.bsmb, original_length, 'nonCG')
                        strand_p_CG = rm_single_CG.record_meth()
                        strand_p_nonCG = rm_single_nonCG.record_meth()
                        loc_dict = LI.Loc_single(read, loc_dict, self.bsmb)
                else:
                    for read in read_inf:
                        rm_pair_CG = CI.PairMethReader(read, strand_p_CG, ref, self.bsmb, original_length, 'CG')
                        rm_pair_nonCG = CI.PairMethReader(read, strand_p_nonCG, ref, self.bsmb, original_length, 'nonCG')
                        strand_p_CG = rm_pair_CG.record_meth()
                        strand_p_nonCG = rm_pair_nonCG.record_meth()
                        loc_dict = LI.Loc_paired(read, loc_dict, self.bsmb)

        if len(self.trim_file) != 0:
            return loc_dict
        else:
            strand_p = [strand_p_CG, strand_p_nonCG]
            return strand_p, loc_dict


    def generator(self):
        '''
        The main method in the class
        1. Draw Mbias plot and generate Mbias table fore per read length in per strand
        2. Decide the trimming positions based on the Mbias plot and generate the trimming file
        3. Show the duplicate reads distribution
        '''

        ref = GR.get_ref(self.ref_file)
        # check: are the input SAM files paired-end or single-end
        strand_p = {}
        if self.single_on:
            strand_p['++'] = {}
            strand_p['-+'] = {}
        else:
            strand_p['++'] = {}
            strand_p['-+'] = {}
            strand_p['+-'] = {}
            strand_p['--'] = {}

        if len(self.trim_file) != 0:
            info("Used the trimming file from the user defined!!")
            info("Ignore both Mbias assessment and trimming decision.")
            loc_dict = self.parser_sambam(strand_p, ref)
            strand_t = self.user_defined_trimming()
        else:
            #modify in 2013-06-04
            strand_t_raw = []
            name_context = [self.name + '_CG', self.name + '_nonCG']
            strand_p, loc_dict = self.parser_sambam(strand_p, ref)
            for i in range(len(strand_p)):
                strand_t_each = self.decide_trim_bp(strand_p[i])
                strand_t_raw.append(strand_t_each)
                MR.mbias_generator(strand_p[i], strand_t_each, name_context[i])
            strand_t = self.decide_final_trimming(strand_t_raw)
            self.produce_final_trim_file(strand_t)
        max_cov = DR.duplicate_report(loc_dict, self.gsize, self.p_poisson, self.name)
        return strand_t, loc_dict, max_cov




