---
title: "Genome_preparation, reads filtering and mapping RNA-Seq from E93-depleted and control embryos samples from *Blatella germanica*"
---

# Genome_preparation, reads filtering and mapping RNA-Seq from E93-depleted and control embryos samples from Blatella germanica
 - Genome and annotation (source Guillem Y):

 -genome: genome/Bgermanica.scaffolds_oneline.fa
 -annotation: genome/BGER_v0.6.2.3.1.gff

 ## 1. Indexing genome:

 - STAR version=2.7.9a

```
STAR --limitGenomeGenerateRAM 40000000000 --runMode genomeGenerate --runThreadN 8 --genomeFastaFiles genome/Bgermanica.scaffolds_oneline.fa --genomeDir genome/genome_indexed_star
```

## 2. Reads preparation
### Data indexes


| sample | index | condition |
|--------|-------|-----------|
| 1 | 16 | Control |
| 2 | 18 | Control |
| 3 | 20 | Control |
| 6 | 23 | dsE93 |
| 8 | 25 | dsE93 |
| 9 | 27 | dsE93 |
| C9 |	15 | Control |
| C11 |19	| Control |
| C12 |	27 | Control |
| T6 |	3	| dsE93 |
| T7 |	5	| dsE93 |
| T14 |	9	| dsE93 |





- md5sum

```
cd results_1 
md5sum *.fastq.gz > md5sum.chk
```

### Join FastQs from same sample

```
cd ..
mkdir joined


cat results_1/1_*_R1_001.fastq.gz > joined/1_R1_001.fastq.gz
cat results_1/1_*_R2_001.fastq.gz > joined/1_R2_001.fastq.gz

cat fastq/2_*_R1_001.fastq.gz > joined/2_R1_001.fastq.gz
cat fastq/2_*_R2_001.fastq.gz > joined/2_R2_001.fastq.gz

cat results_1/3_*_R1_001.fastq.gz > joined/3_R1_001.fastq.gz
cat results_1/3_*_R2_001.fastq.gz > joined/3_R2_001.fastq.gz

cat results_1/6_*_R1_001.fastq.gz > joined/6_R1_001.fastq.gz
cat results_1/6_*_R2_001.fastq.gz > joined/6_R2_001.fastq.gz

cat results_1/8_*_R1_001.fastq.gz > joined/8_R1_001.fastq.gz
cat results_1/8_*_R2_001.fastq.gz > joined/8_R2_001.fastq.gz

cat results_1/9_*_R1_001.fastq.gz > joined/9_R1_001.fastq.gz
cat results_1/9_*_R2_001.fastq.gz > joined/9_R2_001.fastq.gz

cat results_2/C9_*_R1_001.fastq.gz > joined/C9_R1_001.fastq.gz
cat results_2/C9_*_R2_001.fastq.gz > joined/C9_R2_001.fastq.gz

cat results_2/C11_*_R1_001.fastq.gz > joined/C11_R1_001.fastq.gz
cat results_2/C11_*_R2_001.fastq.gz > joined/C11_R2_001.fastq.gz

cat results_2/C12_*_R1_001.fastq.gz > joined/C12_R1_001.fastq.gz
cat results_2/C12_*_R2_001.fastq.gz > joined/C12_R2_001.fastq.gz

cat results_2/T6_*_R1_001.fastq.gz > joined/T6_R1_001.fastq.gz
cat results_2/T6_*_R2_001.fastq.gz > joined/T6_R2_001.fastq.gz

cat results_2/T7_*_R1_001.fastq.gz > joined/T7_R1_001.fastq.gz
cat results_2/T7_*_R2_001.fastq.gz > joined/T7_R2_001.fastq.gz

cat results_2/T14_*_R1_001.fastq.gz > joined/T14_R1_001.fastq.gz
cat results_2/T14_*_R2_001.fastq.gz > joined/T14_R2_001.fastq.gz
```

### Check reads qality

- FastQC v0.11.9
- MultiQC v1.11

```
cd ..
for file in $(ls fastq/joined/*.fastq.gz | grep -v -i "Undetermined")
do
    SAMPLE=`basename $file`
    echo $file
    echo $SAMPLE
    fastqc ${file} -o fastq/Fastqc_reports
done


multiqc  fastq/Fastqc_reports --filename multiqc_report_joinedlanes
```
###  Trimming

- trimmomatic -version 0.39

samples.txt


| sample | index | condition |
|----------|-------------|------|
| 1 | 16 | Control |
| 2 | 18 | Control |
| 3 | 20 | Control |
| 6 | 23 | dsE93 |
| 8 | 25 | dsE93 |
| 9 | 27 | dsE93 |
| C9 |	15 | Control |
| C11 |19	| Control |
| C12 |	27 | Control |
| T6 |	3	| dsE93 |
| T7 |	5	| dsE93 |
| T14 |	9	| dsE93 |


```
python3 trimmomatic.py
```

```
#!/usr/bin/env python

import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threads", metavar="N", help="Number of threads to use", default=12, action="store", type=int, dest="threads")
parser.add_argument("-a", "--adapters", help="Adapter sequences for trimming", default=None, action="store", type=str, dest="adapters")
args = parser.parse_args()


for line in open("samples.txt"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]

    adapters = ""
    if args.adapters != None:    
        adapters = "ILLUMINACLIP:/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:8"


    os.system("trimmomatic PE -threads %d  fastq/joined/%s_R1_001.fastq.gz fastq/joined/%s_R2_001.fastq.gz fastq/cleaned/%s_R1_clened.fastq.gz fastq/cleaned/%s_R1_un.fastq.gz fastq/cleaned/%s_R2_clened.fastq.gz fastq/cleaned/%s_R2_un.fastq.gz %s SLIDINGWINDOW:4:15 2> status/%s.trimming.txt" % (args.threads, sample, sample, sample,sample, sample, sample, adapters, sample))
```

### Check reads qality afrer trimming
```
for file in $(ls fastq/cleaned/*_clened.fastq.gz | grep -v -i "Undetermined")
do
    SAMPLE=`basename $file`
    echo $file
    echo $SAMPLE
    fastqc ${file} -o fastq/cleaned/FastQC
done

multiqc  fastq/cleaned/FastQC --filename multiqc_report_cleaned
```
### remove all unpaird
```
rm fastq/cleaned/*_un.fastq.gz
```

## Mapping reads with STAR

```
mkdir sam_bam

indexname="genome/genome_indexed_star"
outputpath="sam_bam"
fastqdir="fastq/cleaned"

for R1 in $(find $fastqdir -type f  -name "*_R1_clened.fastq.gz"| grep -v "Undetermined")
do

    R2=`echo $R1 | sed 's/_R1_/_R2_/'`
    prebasename=`echo $R1 | sed 's/_R1_clened.fastq.gz//'`
    fname=$(basename $prebasename)

    printf "\n"
        echo Mapping:
	echo $R1
	echo "and"
	echo $R2   
        echo "Output base name:" $fname"_cleaned"	
    printf "\n"

    STAR  --runThreadN 10 \
          --genomeDir $indexname \
          --readFilesIn $R1 $R2 \
          --readFilesCommand zcat\
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix $outputpath/$fname"_clean"


    printf "\n"
    echo  Done with $fname 
    printf "\n"
done;
```
### Mapping Summary

- multiqc 

```
multiqc sam_bam/clean --filename multiqc_report_mapping
```
## Table of counts, batch removal, DEA and enrichment GSEA done in R file:   E93_R_public.Rmd

- R version 4.1.2

- table of counts -in [R] with featureCounts (Rsubread - 2.8.1)
- batch removal -in [R] with SVASeq (sva - 3.42.0) and limma (3.50.0 )
- DEA -in [R] with DeSeq2 (1.34.0)
- GSEA - in [R] clusterProfiler (4.2.1)

- STAR - version 2.7.9a
- trimmomatic -version 0.39
- FastQC v0.11.9
- MultiQC v1.11

 
 
