# RNA_seq
Differential gene expression analysis done on the resistant and susceptible populations of poa annua

# Differential gene expression analysis pipeline 
Differential gene expression is a technique used to find the up- and down-regulated genes in their expression in an experimental design.
Here, we have 6 populations in which 3 are resistant and 3 are susceptible to the herbicide indaziflam

# Tools used 
1.fastp - for removing adapter sequences and quality check of raw fastq files
Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

2.Hisat2 - for aligning the raw reads to the genome
(Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019)

3.SAMTools - for converting and handling SAM/BAM files
(Heng Li and others, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352)

4.EdgeR - for making comparisons between different populations and getting differentially expressed genes.
http://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Procedure

# 1. - Trimming and filtering raw reads
fastp is a command-line tool that can trim the adapter sequences, filter low quality reads and provide quality control files for the raw reads. it provides the quality control files in html and json format, to see the quality of the sequences.

## To install fastp

```
# Create a conda environment 
conda create --name fastp

# installing fastp 
conda install -c bioconda fastp
```

## running fastp

to run fastp, use following command

```
fastp -i ${file}_R1.fastq.gz  -I ${file}_R2.fastq.gz -o ${data}/cleaned_files/${file}_R1.fastq -O ${data}/cleaned_files/${file}_R2.fastq -w 64 \
   --dedup --failed_out ${data}/cleaned_files/fail.fq -j ${data}/cleaned_files/qc/${file}_fastp.json -h ${data}/cleaned_files/qc/${file}_fastp.html 
```
### Meaning of the options used

-i = input of forward strand \
-I = input of reverse strand \
-o = name of forward strand output file, after filtering  \
-O = name of reverse strand output file, after filtering \
-w = number of cores to be used (max it can use = 16) \
--dedup = drop duplicated sequences \
--failed_out = file where all the low quality reads are stored \
-j = file name of quality control file in json format \
-h = file name of quality control file in html format \

for more information, see [fastp github page](https://github.com/OpenGene/fastp)








