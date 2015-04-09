# Introduction #

BSeQC is a quality control package specially for bisulfite sequencing experiments. It can comprehensively evaluate the quality of BS-seq experiments and automatically trim nucleotides with potential technical biases. In addition, BSeQC also support removing duplicate reads and keeping one copy of the overlapping segment in paired-end sequencing.

## BSeQC has the following features: ##


support SAM/BAM files with standard alignment FLAG from all kinds of mapping tools

support WGBS and RRBS

support single-end and paired-end mapping results

support automatical and user-defined trimming

check the biases in diffenrent sequencing length and each mapping strand ('++','-+' for single-end; '++','-+','+-','--' for paired-end)

check CpG and non-CpG

support removing duplicate reads

support keeping one copy of the overlapping segment in paired-end sequencing

support the “[Bis-SNP](http://epigenome.usc.edu/publicationdata/bissnp2011/) strategy” as an alternative method for trimming 5’ bisulfite conversion failure

# Test Datasets #

| **File name** | **Size** | **Discreption** | **Reference** | **QC Figure&table** |
|:--------------|:---------|:----------------|:--------------|:--------------------|
|[mNPC\_Paired\_BSseq\_chr1](http://compbio.tongji.edu.cn/~linxq/mNPC_Paired_BSseq_chr1.tar.gz) | 2.2G | mouse ES-derived neuronal progenitor cells (GSE30202)| Stadler, M.B., et al. (2011)|[mNPC\_Paired\_BSseq\_chr1\_BSeQC\_report](http://compbio.tongji.edu.cn/~linxq/mNPC_Paired_BSseq_chr1_BSeQC_figure_table.tar.gz)|
|  |  |  |  | [mNPC\_Paired\_BSseq\_chr1\_Bis-SNP\_strategy\_report](http://compbio.tongji.edu.cn/~linxq/mNPC_Paired_BSseq_chr1_BSeQC_Bis-SNP_strategy_result.tar.gz)|
|[H1\_Single\_BSseq\_chr1\_replicate1](http://compbio.tongji.edu.cn/~linxq/H1_Single_Methylseq_chr1_replicate1.tar.gz) | 4.9G | H1 cell line (GSE16256) | Lister, R., et al. (2009)| [H1\_Single\_BSseq\_chr1\_replicate1\_BSeQC\_report](http://compbio.tongji.edu.cn/~linxq/H1_Single_Methylseq_chr1_replicate1_BSeQC_figure_table.tar.gz)|
|[H1\_Single\_BSseq\_chr1\_replicate2 ](http://compbio.tongji.edu.cn/~linxq/H1_Single_Methylseq_chr1_replicate2.tar.gz)| 5.9G | H1 cell line (GSE16256) | Lister, R., et al. (2009)| [H1\_Single\_BSseq\_chr1\_replicate2\_BSeQC\_report](http://compbio.tongji.edu.cn/~linxq/H1_Single_Methylseq_chr1_replicate2_BSeQC_figure_table.tar.gz)|



