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
-h = file name of quality control file in html format 

for more information, see [fastp github page](https://github.com/OpenGene/fastp)

# 2 Aligning the raw reads 
After initial filtering, trimming, and quality control of the raw reads, the next step is to align the raw reads. We used HISAT2 to align the raw reads to the reference transcriptome.

## Building the genome index
Before running the alignment, hisat needs trancriptome to be indexed. 

```
# Building the trancriptome index 
hisat2-build -p 64 data/run/maheym/poa_annua/poa_annua_cds_index/PoaAn.maker.cds \
/data/run/maheym/poa_annua/poa_annua_cds_index/PoaAn.maker.index 
```
where
-p = number of cores to use \
data/run/maheym/poa_annua/poa_annua_cds_index/PoaAn.maker.cds = location of the trancriptome sequence \
/data/run/maheym/poa_annua/poa_annua_cds_index/PoaAn.maker.index = location and name of indexed transcriptome to be used for alignment

## running the Hisat2 aligner
To run the aligner, we need
a. cleaned forward reads (made in previous using fastp)
b. cleaned reverse reads (made in previous using fastp)
c. Indexed reference transcriptome to be aligned (if reference transcriptome is unavailable, consider using Trinity for de-novo assemby or Salmon, kallisto for transcript quantification without aligning)

```
hisat2 -p 70 --quiet \
               -x /data/run/maheym/poa_annua/poa_annua_cds_index/PoaAn.maker.index \
               -1 /data/run/maheym/poa_annua/raw_data/${file}_1.fq.gz \
               -2 /data/run/maheym/poa_annua/raw_data/${file}_2.fq.gz \
               -S /data/run/maheym/poa_annua/Poa.maker.transcripts_${file}.sam
```
where
-p = number of cores to use
--quiet = to prevent priting to terminal, except sequences or serious errors
-x = input indexed reference transcriptome
-1 = forward reads 
-2 = reverse reads
-S = output as SAM format (default output is stdout)

------------------------
### NOTE: This part is specific to the server we used and has no relation with aligning or the rna-seq pipeline. Everything can be run without this part of code.
The script was run on the server using torque. #PBS are various flags that are needed by job batcher. 
```
#PBS -N HiSat2_transcripts_poa_annua
#PBS -l nodes=1:ppn=70,mem=501gb
#PBS -d /data/run/maheym/poa_annua
#PBS -M maheymoh@msu.edu
#PBS -m abe

# Loading the modules
module load HISAT2
module load SAMTools
```
-N = name of the job \
-l = numer of nodes, cpus and memory required \
-d = working directory \
-M = email id to get updates
-m = abe means update, if the job is "a"borts,b"egins,"e"nd 

installing HISAT2
```
conda install -c bioconda hisat2
```
--------------------------

# 3 Sorting and indexing aligned SAM outputs
After aligning with HISAT2, the reads are un-sorted. The reads needs to be sorted for downward analysis. after sorting the reads needs to be indexed. The indexing is process of getting the position of the reads in alignment file. Due to large size of the files, indexing helps the softwares to run faster and efficiently. SAM files can go upto 50GB per file, depending on size of your transcriptome. Thus it is better to convert SAM (sequence alignment mapping) to BAM (binary alignement mapping) to save space and run downward analysis faster.

## Converting SAM to BAM and sorting the BAM files

```
samtools view -hb -@ 35 Poa.maker.transcripts_${file}.sam | samtools sort -@ 35 \
         -o /data/run/maheym/poa_annua/poa_transcripts/Poa.maker.transcripts_${file}.bam
```

view = flag to convert SAM to BAM \
-hb(to flags combined in one, similar to -h, -b) = -h = includes headers, -b - output as BAM \
| = piping the ouput \
-@ = number of cores to use \
-sort = sort the reads in BAM file \ 
-o = output location for sorted BAM file

## Indexing the BAM files
After sorting we can index the BAM files. 
```
samtools index /data/run/maheym/poa_annua/poa_transcripts/Poa.maker.transcripts_${file}.bam
        rm Poa.maker.transcripts_${file}.sam
```
after converting SAM to BAM, we can remove SAM files to save space. We can always convert BAM to SAM or vice-versa. They are practically the same information, just SAM is in human readable form (text)  while BAM is in binary form. 

# getting raw reads count 












