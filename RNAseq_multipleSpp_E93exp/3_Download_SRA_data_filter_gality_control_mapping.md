- Author: Gabriela Machaj
- Year: 2021

# RNA-seq data analysis

This script describes: RNA-seq data retrieval, the QC, reads trimming, and read mapping.


## RNA-seq data downloading with fastq-dump

- List of species, and selected SRR samples in **species_list** .


```
python  prefetch.py 
```

```
#!/usr/bin/env python
import os
import sys



for line in open("species_list"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]

    os.system("prefetch -v %s" % (sample))
    os.system("fastq-dump --outdir fastq --split-files sra/%s.sra" % (sample))
```

## Check quality of raw fastq files

```
for file in $(ls fastq/*.fastq | grep -v -i "Undetermined")
do
    SAMPLE=`basename $file`
    echo $file
    echo $SAMPLE
    fastqc ${file} -o fastq/FastQC_raw
done
```

### Summarize QC

```
multiqc FastQC_raw --interactive --force --filename multiqc_raw
```

## Trim reads with trim_galore

- for single end reads

```
python trimgalore_single.py
```

```
#!/usr/bin/env python
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threads", metavar="N", help="Number of threads to use", default=8, action="store", type=int, dest="threads")

args = parser.parse_args()


for line in open("species_list"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]

    os.system("TrimGalore-0.6.6/trim_galore  --phred33 --trim-n --cores %d --output_dir fastq/trim_galore fastq/%s_1.fastq" % (args.threads, sample) )
```

- for paired end reads

```
python trimgalore_paired.py
```

```
#!/usr/bin/env python
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threads", metavar="N", help="Number of threads to use", default=8, action="store", type=int, dest="threads")

args = parser.parse_args()


for line in open("SRA_list"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]

    os.system("TrimGalore-0.6.6/trim_galore --paired --phred33 --trim-n --cores %d  --output_dir fastq/trim_galore fastq/%s_1.fastq.gz fastq/%s_2.fastq.gz" % (args.threads, sample, sample) )
```

### Check quality after QC

```
for file in $(ls trim_galore/*.fq | grep -v -i "Undetermined")
do
    SAMPLE=`basename $file`
    echo $file
    echo $SAMPLE
    fastqc ${file} -o trim_galore/FastQC_tg
done
```

### Summarize QC

```
multiqc trim_galore/FastQC_tg --interactive --force --filename multiqc_trim_galore
```

## Reads mapping

- for paired end reads
```
python mapping_paired.py
```

```
#!/usr/bin/env python
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threads", metavar="N", help="Number of threads to use", default=12, action="store", type=int, dest="threads")

args = parser.parse_args()

for line in open("SRA_list"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]


    read_files = None
        read_files = "SRA/fastq/trim_galore/%s_1_val_1.fq SRA/fastq/trim_galore/%s_2_val_2.fq" % (sample, sample)

    os.system("STAR --runThreadN %d --genomeDir data/species_name/genome_indexed_star --readFilesIn %s --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/%s/" % (args.threads, read_files, sample) )

```

- for single end reads
```
python mapping_single.py
```

```
#!/usr/bin/env python
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threads", metavar="N", help="Number of threads to use", default=12, action="store", type=int, dest="threads")
args = parser.parse_args()


for line in open("SRA_list"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]


    read_files = None
        read_files = "SRA/fastq/trim_galore/%s_1_trimmed.fq" % (sample,)


    os.system("STAR --runThreadN %d --genomeDir data/species_name/genome_indexed_star --readFilesIn %s --outSAMtype BAM SortedByCoordinate --outFileNamePrefix bam/%s/" % (args.threads, read_files, sample) )
```

### Summarize mapping

```
multiqc bam --interactive --force --filename multiqc_mapping
```
## Programs and versions:

- STAR - v.2.7.9a
- trim-galore v.0.6.7 
- MultiQC v.1.11
- FastQC v.0.11.9


### Next step 

- Count reads mapped to genes --> [4_R_featurecount.Rmd](https://github.com/ylla-lab/Embryonic_E93/blob/master/RNAseq_multipleSpp_E93exp/4_R_featurecount.Rmd)
