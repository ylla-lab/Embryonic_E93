---
title: "Genome_preparation-INSECTS-E93"
author: "Gabriela Machaj"
date: "11/16/2021"
output: html_document
---

- Author: Gabriela Machaj
- Date:  09/03/2022

## Data dowloadadet frm NCBI (https://www.ncbi.nlm.nih.gov/genome/?term=) 16.11.2021r 

- species directory list:

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
------
 - B.germanica - source G.Y.
------
## Data dowloadadet frm NCBI (https://www.ncbi.nlm.nih.gov/genome/?term=) 31.01.2022r

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
_____________________________________

 - G.bimaculatus   - suurce G.Y.

_________________________________

- in each directory data renemed as:

  - genome: genome.fna
  - annotation: annotation.gff

### Index genome

- for each species except B.germanica and D.citri (%s = species_directory_name) for B.germanica and D.citri (B_D_names)

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
- for B.grmanica and D.citi genomes

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
