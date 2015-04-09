# Introduction #

Bisulfite sequencing (BS-seq) has been extensively used for DNA methylome study. Ideally, BS-seq experiment should be able to directly and exactly identify the methylation state of a DNA fragment from the original genome. However, current BS-seq protocols still possess several intrinsic biases, which will impact methylation level estimation, such as overhang end-repair, 5’ bisulfite conversion failure, sequencing into the adaptor and 3’ low sequencing quality.  Since BS-seq experiments are widely used and the resulting data will continue to grow exponentially in the near future, there is a strong need for a dedicated QC tool to evaluate and remove potential technical biases in BS-seq experiments. Here, we developed BSeQC package. It can comprehensively evaluate the quality of BS-seq experiments and automatically trim nucleotides with potential technical biases that may result in inaccurate methylation estimation. In addition, BSeQC also support removing duplicate reads and keeping one copy of the overlapping segment in paired-end sequencing.



# Installation #

## Prerequisite: ##

  * Python>=2.6
  * Numpy
  * Scipy
  * Matplotlib  (We recommend to install this python package. If not, we only can produce the Mbias table, not the Mbias plot)

## Install BSeQC in defaut location ##
```
1. tar zxf BSeQC-VERSION.tar.gz
2. cd BSeQC-VERSION
3. python setup.py install  
 #Skip step4 if '/usr/local/lib/python2.7/site-packages' is already included in your PYTHONPATH.
4. export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH 
 #Skip step5 if '/usr/local/bin' is already included in your PATH.
5. export PATH=/user/local/bin:$PATH
 #To make permanent change to your PYTHONPATH or PATH variables, copy the commands (step4 and step5) into your
 #'/home/user/.bashrc' or '/home/user/.bash_profile'. 
```
## Install BSeQC in user specified  location ##
```
1. tar zxf BSeQC-VERSION.tar.gz
2. cd BSeQC-VERSION
 #You need to change '/home/user/' accordingly
3. python setup.py install --prefix=/home/user/
 #setup PYTHONPATH, so that BSeQC knows where to import modules
4. export PYTHONPATH=/home/user/lib/pythonX.Y/site-packages:$PYTHONPATH
 #setup PATH, so that system knows where to find executable files. 
5. export PATH=/home/user/bin:$PATH 
 #To make permanent change to your PYTHONPATH or PATH variables, copy the commands (step4 and step5) into your
 #'/home/user/.bashrc' or '/home/user/.bash_profile'. 
```


# Usage Information #

mbias: BSeQC main executable function

| **Options** | **Description** | **Requirement** |
|:------------|:----------------|:----------------|
| -s or  --sam | The SAM file for quality analysis; Multiple SAM file should be separated by the ','. | Mandatory |
| -r or  --ref | The reference genome fasta file. | Mandatory |
| -t or  --samtools | The path of samtools. | Optional |
| -n or  --name | The name for the output plot and table. Default = 'NA'. | Optional |
| -l or  --len	| If the original mapping reads have been trimmed with adapter or other reasons, the original read length for the sam file should be set. Multiple length can also be separated by ','. If the read length of two mates in paired-end is different, please separated by '-'. | Optional |
| -p or --pvalue | How many stds will be set for the trimming cutoff. Default = 2. |	 Optional |
| --drift | How many drifts(%mC) will be set for the trimming cutoff. Default = 2. | Optional |
| -f or --trim\_file | User can determine the trimming bp by the trim file. | Optional |
| -a or --auto	| Automatically trim the biased bp. If not you can use the Mcall biases plot to manually decide how many bps to trim and make a trimming file. Default = True. | Optional |
|-o or --remove\_overlap| Keep only one copy of the overlapping segment of two read mates in paired-end seq. Default = True.  |Optional|
|--filter\_dup|Remove duplicate reads resulting from possible over-amplification. Default = True. |Optional|
|--p\_poisson| Pvalue cutoff Poisson distribution test in removing duplicate reads. Default = 1e-5.| Optional |
|-g or --gsize | Effective genome size for calculate max coverage. It can be 1.0e7 or 10000000, or shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8). It is restricted by --filter\_dup.| Optional |
|--not\_mapping|Whether keep the not-unique mapping reads in the filter SAM file. Default = True. | Optional |


nonuniform: use the “Bis-SNP strategy” as an alternative method for trimming 5’ bisulfite conversion failure
| **Options** | **Description** | **Requirement** |
|:------------|:----------------|:----------------|
| -s or  --sam | The SAM file for quality analysis; Multiple SAM file should be separated by the ','. | Mandatory |
| -r or  --ref | The reference genome fasta file. | Mandatory |
| -t or  --samtools | The path of samtools. | Optional |
| -n or  --name | The name for the output plot and table. Default = 'NA'. | Optional |
|--filter\_dup|Remove duplicate reads resulting from possible over-amplification. Default = True. |Optional|
|-o or --remove\_overlap| Keep only one copy of the overlapping segment of two read mates in paired-end seq. Default = True. |Optional|
|--p\_poisson| Pvalue cutoff Poisson distribution test in removing duplicate reads.| Optional |
|-g or --gsize | Effective genome size for calculate max coverage. It can be 1.0e7 or 10000000, or shortcuts: 'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8). It is restricted by --filter\_dup.| Optional |
|--not\_mapping|Whether keep the not-unique mapping reads in the filter SAM file. Default = True. | Optional |



### Strand symbol and information: ###
|Read type|	Strand Symbol|	Strand information|
|:--------|:-------------|:------------------|
|Single-end|	++|	Watson strand|
|Single-end|	-+|	Crick strand|
|Paired-end|	++|	Forward strand of Watson of reference (mate1)|
|Paired-end|	+-|	Reverse strand of Watson of reference (mate2)|
|Paired-end|-+|	Forward strand of Crick of reference (mate1)|
|Paired-end|--|	Reverse strand of Crick of reference (mate2)|


# Example #

**mbias: BSeQC main executable function**
```
##paired-end (test data: mNPC_chr1)
#just focus on the mbias trimming
bseqc mbias -s SRR299067_chr1.bam,SRR299069_chr1.bam,SRR299071_chr1.bam,SRR299073_chr1.bam,SRR299075_chr1.bam,SRR299077_chr1.bam,SRR299079_chr1.bam,SRR299081_chr1.bam,SRR299083_chr1.bam,SRR299068_chr1.bam,SRR299070_chr1.bam,SRR299072_chr1.bam,SRR299074_chr1.bam,SRR299076_chr1.bam,SRR299078_chr1.bam,SRR299080_chr1.bam,SRR299082_chr1.bam -r mm9.fa -n mNPC_Paired_BSseq_chr1 -o --filter_dup

#filter duplicate reads, remove one copy of the overlapping segment, and not keep not_unique mapping reads
bseqc mbias -s SRR299067_chr1.bam,SRR299069_chr1.bam,SRR299071_chr1.bam,SRR299073_chr1.bam,SRR299075_chr1.bam,SRR299077_chr1.bam,SRR299079_chr1.bam,SRR299081_chr1.bam,SRR299083_chr1.bam,SRR299068_chr1.bam,SRR299070_chr1.bam,SRR299072_chr1.bam,SRR299074_chr1.bam,SRR299076_chr1.bam,SRR299078_chr1.bam,SRR299080_chr1.bam,SRR299082_chr1.bam -r mm9.fa -n mNPC_Paired_BSseq_chr1_roverlap_dup -g 1.9e8 --not_mapping

#use trimming file 
bseqc mbias -s SRR299067_chr1.bam,SRR299069_chr1.bam,SRR299071_chr1.bam,SRR299073_chr1.bam,SRR299075_chr1.bam,SRR299077_chr1.bam,SRR299079_chr1.bam,SRR299081_chr1.bam,SRR299083_chr1.bam,SRR299068_chr1.bam,SRR299070_chr1.bam,SRR299072_chr1.bam,SRR299074_chr1.bam,SRR299076_chr1.bam,SRR299078_chr1.bam,SRR299080_chr1.bam,SRR299082_chr1.bam -r mm9.fa -n mNPC_Paired_BSseq_chr1_roverlap_dup -g 1.9e8 --not_mapping -f mNPC_Paired_BSseq_chr1_trim_file.txt


##single-end
#just focus on on the mbias trimming (testdata: H1_chr1)
#for replicate1: set the samtools path; set the original sequence length, because some bps in the 3' end of the reads with low quality or adapter have been trimmed during mapping   
bseqc mbias -s methylC-seq_H1_r1_noaq_chr1_43.bam,methylC-seq_H1_r1_noaq_chr1_52.bam,methylC-seq_H1_r1_noaq_chr1_53.bam,methylC-seq_H1_r1_noaq_chr1_76.bam,methylC-seq_H1_r1_noaq_chr1_87.bam,methylC-seq_H1_r1_noaq_chr1_88.bam -r hg19.fa -l 43,52,53,76,87,88 -n H1_Single_BSseq_chr1_replicate1 --filter_dup

#for replicate2
bseqc mbias -s methylC-seq_H1_r2_noaq_chr1.bam -r hg19.fa  -n H1_Single_BSseq_chr1_replicate2  -l 87 --filter_dup

```

**nonuniform: use the “Bis-SNP strategy” as an alternative method for trimming 5’ bisulfite conversion failure**
```
##test data: mNPC_chr1
bseqc nonuniform -s SRR299067_chr1.bam,SRR299069_chr1.bam,SRR299071_chr1.bam,SRR299073_chr1.bam,SRR299075_chr1.bam,SRR299077_chr1.bam,SRR299079_chr1.bam,SRR299081_chr1.bam,SRR299083_chr1.bam,SRR299068_chr1.bam,SRR299070_chr1.bam,SRR299072_chr1.bam,SRR299074_chr1.bam,SRR299076_chr1.bam,SRR299078_chr1.bam,SRR299080_chr1.bam,SRR299082_chr1.bam -r /mnt/Storage/data/Bowtie/mm9.fa -n Bis-SNP_strategy --filter_dup -o --not_mapping
```

**rrbs: support RRBS data**
```
##RRBS paired-end 
#test data: SRR726536
bseqc rrbs -s SRR726536_chr1.sam -r mm9.fa -n rrbs_p
##RRBS single-end
#test data: SRR788619 & SRR788620
bseqc rrbs -s SRR788619.sam -r mm9.fa -n rrbs_s_rep1
```

# Output #

| **name** | **content** |
|:---------|:------------|
|**mbias: BSeQC main executable function**|  |
|name + 'CG\_Mbias\_plot.pdf'| the CpG Mbias plot for each read length in each stand|
|name + 'nonCG\_Mbias\_plot.pdf'| the non-CpG Mbias plot for each read length in each stand|
|name + '_Dup\_dis.pdf'_|    the duplicate reads distribution (when --filter\_dup be set true) |
|name\_CG\_strand\_readlength.csv(in directory:name+'_Mbias\_table)_| the CpG Mbias table for each read length in each stand|
|name\_nonCG\_strand\_readlength.csv(in directory:name+'_Mbias\_table)_| the non-CpG Mbias table for each read length in each stand|
|name\_CG\_trim\_file.txt(in directory:name+_Mbias\_table)_| the trimming decision for each read length in each stand from the CpG Mbias plot|
|name\_nonCG\_trim\_file.txt(in directory:name+_Mbias\_table)_| the trimming decision for each read length in each stand from the non-CpG Mbias plot|
|name\_final\_trim\_file.txt | the most stringent trimming decision made by either CpG or non-CpG cytosines M-bias plots|
|name\_BSeQC\_mbias\_filter\_report.txt | the detail trimming and filtering report |
|**nonuniform: use the “Bis-SNP strategy” as an alternative method for trimming 5’ bisulfite conversion failure**|  |
|name\_trimmed\_nucleotides\_dis.pdf | the distribution of the number of trimmed nucleotides based on the “Bis-SNP strategy”|
|name\_BSeQC\_mbias\_filter\_report.txt | the detail trimming and filtering report |

# NOTE #
Because different DNA strands and read lengths can have distinct biases, BSeQC will trim them differently. If users are concerned about the different coverage on two strands, we can use the parameter  (-f TRIM\_FILE) to let Watson(+) and Crick(-) strands be trimmed with same nucleotides in the 5’ and 3’.

### The format of the trim\_file(split by tab): ###
|strand|	read\_length |	trim\_5'_pos_|	trim\_3'_pos_|
|:-----|:-------------|:--------|:--------|
|++ |	100 |	1|	98|
|++|	88|	2|	78|
|-+|	100|	2|	98|
|-+|	87|	2|	78|