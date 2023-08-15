# RNA_seq
Differential gene expression analysis done on the resistant and susceptible populations of poa annua

# Differential gene expression analysis pipeline 
Differential gene expression is a technique used to find the up- and down-regulated genes in their expression in an experimental design.
Here, we have 6 populations in which 3 are resistant and 3 are susceptible to the herbicide indaziflam

# tools used 
1. SAMTools - for converting and handling SAM/BAM files
(Heng Li and others, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352)

2. Hisat2 - for aligning the raw reads to the genome
(Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019)

3. fastp - for removing adapter sequences and quality check of raw fastq files
Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560


