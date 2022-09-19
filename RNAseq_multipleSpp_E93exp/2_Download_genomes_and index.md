
# Download and Index genomes

- Author: Gabriela Machaj
- Date:  11/16/2021

## Data downloaded from NCBI - 16.11.2021

https://www.ncbi.nlm.nih.gov/genome/?term= + SPECIES

Species directory list:

- A. mellifera 
- A. pisum
- A. stephensi
- B. mori
- C. dipterum
 - C. floridanus
 - C. quinquefasciatus
 - D. melanogaster
 - H. saltator
 - M. sexta
 - N. ugens
 - T. castaneum
 - T. sarcophagae
 - A. gossypii
 - B. dorsalis
 - B. oleae
 - B. tabaci
 - D. citri
 - F. candida
 - H. illucens
 - L. heterotoma
 - M. genalis
 - M. pharaonis
 - O. brunneus
 - P. polytes
 - P. rapae
 - P. xylostella
 - S. exigua
 - S. frugiperda
 - V. tameamea
 - Z. cucurbitae
 - Z. nevadensis
_____________________________________

Genomes not in NCBI:

 - G. bimaculatus   - source G.Y.
 - B. germanica   - source G.Y.

_________________________________

- in each directory data renemed as:

  - genome: genome.fna
  - annotation: annotation.gff

## Index genomes

- for each species (except *B. germanica* and *D. citri*) (%s = species_directory_name). 

```
#!/usr/bin/env python
import os
import sys

for line in open("species_directory_name"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]

    os.system("STAR --limitGenomeGenerateRAM 40000000000 --runMode genomeGenerate --runThreadN 16 --genomeFastaFiles %s/genome.fna --genomeDir %s/genome_indexed_star" % ( sample, sample) )
```
- for *B. germanica* and *D. citri* genomes

```
#!/usr/bin/env python
import os
import sys

for line in open("B_D_names"):
    if line.startswith("#"): continue
    line = line.rstrip("\r\n").split()
    sample = line[0]

    os.system("STAR --limitGenomeGenerateRAM 400000000000 --genomeChrBinNbits 12 --genomeSAindexNbases 12  --runMode genomeGenerate --runThreadN 12 --genomeFastaFiles %s/genome.fna --genomeDir %s/genome_indexed_star"  % ( sample, sample) )
```

- STAR version: 2.7.9a
