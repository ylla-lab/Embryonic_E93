---
title: "Download_SRA_data_filter_gality_control_mapping-INSECTS-E93"
author: "Gabi_M"
output: html_document
---
## Data 
Filtered in R and manually SRA library names were listed  separatelly  for each example species_list 

B.germanica

```
SRR5656360
SRR5656359
SRR5656362
SRR5656368
SRR5656367
SRR5656366
SRR5656363
SRR5656361
SRR5656364
SRR5656365
SRR6784708
SRR6784709
```
#### Species SRA lists:

- A.melliferan 
- A.pisum
- A.stephensi
- B.mori
- C.dipterum
 - C.floridanus
 - C.quinquefasciatus
 - D.melanogaster
 - H.saltator
 - M.sexta
 - N.ugens
 - T.castaneum
 - T.sarcophagae
 - A.gossypii
 - B.dorsalis
 - B.oleae
 - B.tabaci
 - D.citri
 - F.candida
 - H.illucens
 - L.heterotoma
 - M.genalis
 - M.pharaonis
 - O.brunneus
 - P.polytes
 - P.rapae
 - P.xylostella
 - S.exigua
 - S.frugiperda
 - V.tameamea
 - Z.cucurbitae
 - Z.nevadensis
_-  G.bimaculatus

## Data downloading and fastq-dump
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

### Check quality of raw fastq files
```
for file in $(ls fastq/*.fastq | grep -v -i "Undetermined")
do
    SAMPLE=`basename $file`
    echo $file
    echo $SAMPLE
    fastqc ${file} -o fastq/FastQC_raw
done
```

#### Summarize QC
```
multiqc FastQC_raw --interactive --force --filename multiqc_raw
```
### Trim reads with trim_galore

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

#### Summarize QC

```
multiqc trim_galore/FastQC_tg --interactive --force --filename multiqc_trim_galore
```

### Reads mapping

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

#### Summarize mapping

```
multiqc bam --interactive --force --filename multiqc_mapping
```
#### Programs and versions:

- STAR - v.2.7.9a
- trim-galore v.0.6.7 
- MultiQC v.1.11
- FastQC v.0.11.9


### Next step count reads mapped to genes --> 4_R_featurecount.Rmd