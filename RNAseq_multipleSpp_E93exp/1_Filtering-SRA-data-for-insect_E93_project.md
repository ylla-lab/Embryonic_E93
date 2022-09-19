# Project: filtering SRA data for insects to analyze E93 expression in embryo/egg vs. pupa (holometabolous) and embryo/egg vs. nymph/juvenile (a-, neo- and himimetabolous).

Author: Gabriela Machaj
Year: 2021

## List of annotated genomes

List of isects genomes with annotation (filter: annotated) was
downloaded from:
<https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=50557> (accesed:
8.11.2021) (354- total) – \> genome-list 240

-   Filename: annotated_genomes_taxa50557

List of isects genomes with annotation (filter: annotated) was
downloaded from:
<https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=50557> once again
( accesed: 25.01.2022) (384) – \> genome list 260

-   Filename: annotated_genomes_taxa50557_01_2022

List of hexapoda hexapoda genomes with annotation (filter: annotated)
was downloaded from:
<https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=6960> ( accesed:
25.01.2022) (388) –\> genome list 263

-   filename: annotated_genomes_hexapoda

Gryllus bimaculatus - not annotated in NCBI - annotation provided by
G.Ylla was manually added to the list

## Download and filter SRA mtetadata

Metadata for SRA Experiments was first filtered in brwser
(<https://www.ncbi.nlm.nih.gov/sra>) using querys and filters:

### For insects:

1.  Embryo

    -   query: ((txid50557\[Organism:exp\]) AND embryo OR egg) NOT
        mirna-seq
    -   filters: RNA, Illumina, fastq, RNA-SEq (91 data) /Other (8,644
        data ) in RunSelector:

-   filename: embryo_egg_RNA_Seq.txt
-   filename: embryo_egg_Other.txt

1.  Nymph
    -   query: ((txid50557\[organism:exp\]) AND nymph OR juvenile ) NOT
        mirna-seq

    -   filters: RNA, Illumina, fastq, RNA-SEq (48 data) /Other (16 155
        data) in RunSelector:

-   filename: nymph_juvenile_OTHER.txt
-   filename: nymph_juvenile_RNA_Seq.txt

1.  Pupa
    -   query: ((txid50557\[Organism:exp\]) AND pupa) NOT mirna-seq

    -   filters: RNA, Illumina, fastq, RNA-SEq (0 data) /Other (1058
        data) in RunSelector:

-   filename: pupa_OTHER.txt

### For hexapoda:

1.  Embryo
    -   query: ((txid6960\[Organism:exp\]) AND embryo OR egg) NOT
        mirna-seq

    -   filters: RNA, Illumina, fastq, RNA-SEq (91 data) /Other (9,834
        data) in RunSelector

-   filename: embryo_egg_hexa_RNA_Seq.txt
-   filename: embryo_egg_hexa_OTHER.txt

1.  Nymph
    -   query: ((txid6960\[Organism:exp\]) AND nymph OR juvenile ) NOT
        mirna-seq

    -   filters: RNA, Illumina, fastq, RNA-SEq (48 data) /Other (16 155
        data) in RunSelector:

-   filename: nymph_juvenile_hexa_OTHER.txt - filename:
    nymph_juvenile_hexa_RNA_Seq.txt

<!-- -->

    !!!!! there is an error in SRA, if you select "OTHER" there is RNA-Seq, if you select "RNA-Seq" there is "other and RIP -Seq!!!!!!!!!

### prossesing data in R to select all available SRA experiments:

The metadata sets downloaded from the SRA database were first filtered
to leave only those for annotated species, then the metadata for
embryo/egg were compared with the metadata for pupa and nymph/juvenile –
\> we selected all available data for annotated species containing
RNA-Seq data for both developmental stages (embryo/egg + nymph/juvenile
and embryo/egg + pupa)

\### Read data and R (v.4.1.3)

# read new data for juvenile and egg

``` r
SRA_nymph_juvenile_OTHER<-read_csv("nymph_juvenile_OTHER.txt")
nrow(SRA_nymph_juvenile_OTHER)

SRA_nymph_juvenile_RNA_Seq<-read_csv("nymph_juvenile_RNA_Seq.txt")
nrow(SRA_nymph_juvenile_RNA_Seq)


SRA_nymph_juvenile <- bind_rows(SRA_nymph_juvenile_OTHER,SRA_nymph_juvenile_RNA_Seq)
nrow(SRA_nymph_juvenile)

SRA_embryo_egg_OTHER<-read_csv("embryo_egg_OTHER.txt")
nrow(SRA_embryo_egg_OTHER)

SRA_embryo_egg_RNA_Seq<-read_csv("embryo_egg_RNA_Seq.txt")
nrow(SRA_embryo_egg_RNA_Seq)

SRA_nymph_juvenile_hexa_OTHER<-read_csv("nymph_juvenile_hexa_OTHER.txt")

SRA_nymph_juvenile_hexa_RNA_Seq<-read_csv("nymph_juvenile_hexa_RNA_Seq.txt")

SRA_nymph_juvenile_hexa <- bind_rows(SRA_nymph_juvenile_hexa_OTHER,SRA_nymph_juvenile_hexa_RNA_Seq)

SRA_embryo_egg_hexa_OTHER<-read_csv("embryo_egg_hexa_OTHER.txt")

SRA_embryo_egg_hexa_RNA_Seq<-read_csv("embryo_egg_hexa_RNA_Seq.txt")
```

### filter SRA data with annotated genome for new data set - juvenile egg

``` r
genome_list<- genome %>% select(Organism.Name) %>% distinct(Organism.Name) # list of org with annotated genomes
nrow(genome_list)

genome_list_new<- genome_new %>% select(Organism.Name) %>% distinct(Organism.Name) # list of org with annotated genomes
nrow(genome_list_new)

genome_list_hexa<- genome_hexa %>% select(Organism.Name) %>% distinct(Organism.Name) # list of org with annotated genomes
nrow(genome_list_hexa)

embryo_egg_genome_OTHER<-filter(SRA_embryo_egg_OTHER,Organism %in% genome_list_new$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(embryo_egg_genome_OTHER)

embryo_egg_genome_RNA_Seq<-filter(SRA_embryo_egg_RNA_Seq,Organism %in% genome_list_new$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(embryo_egg_genome_RNA_Seq)

embryo_egg_genome_OTHER %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_genome_OTHER_list<-embryo_egg_genome_OTHER %>% select(Organism) %>% distinct(Organism)

embryo_egg_genome_RNA_Seq %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_genome_RNA_Seq_list<-embryo_egg_genome_RNA_Seq %>% select(Organism) %>% distinct(Organism)

nymph_juvenile<-filter(SRA_nymph_juvenile,Organism %in% genome_list_new$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(nymph_juvenile)

nymph_juvenile %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_juvenile_list<-nymph_juvenile %>% select(Organism) %>% distinct(Organism)

# for hexapodes

embryo_egg_genome_hexa_OTHER<-filter(SRA_embryo_egg_hexa_OTHER,Organism %in% genome_list_hexa$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(embryo_egg_genome_hexa_OTHER)

embryo_egg_genome_hexa_RNA_Seq<-filter(SRA_embryo_egg_hexa_RNA_Seq,Organism %in% genome_list_hexa$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(embryo_egg_genome_hexa_RNA_Seq)

embryo_egg_genome_hexa_OTHER %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_genome_hexa_OTHER_list<-embryo_egg_genome_hexa_OTHER %>% select(Organism) %>% distinct(Organism)

embryo_egg_genome_hexa_RNA_Seq %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_genome_hexa_RNA_Seq_list<-embryo_egg_genome_hexa_RNA_Seq %>% select(Organism) %>% distinct(Organism)

nymph_juvenile_hexa<-filter(SRA_nymph_juvenile_hexa,Organism %in% genome_list_hexa$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(nymph_juvenile_hexa)

nymph_juvenile_hexa %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_juvenile_hexa_list<-nymph_juvenile_hexa %>% select(Organism) %>% distinct(Organism)
```

### filter SRA data with annotated genome (firsta time -no hexapoda no egg juvenile)

``` r
genome_list<- genome %>% select(Organism.Name) %>% distinct(Organism.Name) # list of org with annotated genomes
nrow(genome_list)

embryo_genome<-filter(SRA_embryo,Organism %in% genome_list$Organism.Name) # filter Embryo etc. data by Organism  present in genome
nrow(embryo_genome)
embryo_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_genome_list<-embryo_genome %>% select(Organism) %>% distinct(Organism)

pupa_genome<-filter(SRA_pupa,Organism %in% genome_list$Organism.Name) 
nrow(pupa_genome)
pupa_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()
pupa_genome_list<-pupa_genome %>% select(Organism) %>% distinct(Organism)


pupa_genome_hexa<-filter(SRA_pupa,Organism %in% genome_list_hexa$Organism.Name) 
nrow(pupa_genome_hexa)
pupa_genome_list_hexa<-pupa_genome_hexa %>% select(Organism) %>% distinct(Organism)


nymph_genome<-filter(SRA_nymph,Organism %in% genome_list$Organism.Name) 
nrow(nymph_genome)
nymph_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_genome_list<-nymph_genome %>% select(Organism) %>% distinct(Organism)

larva_genome<-filter(SRA_larva,Organism %in% genome_list$Organism.Name) 
nrow(larva_genome)
larva_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

larva_genome_list<-larva_genome %>% select(Organism) %>% distinct(Organism)

# for hexa

embryo_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_genome_list<-embryo_genome %>% select(Organism) %>% distinct(Organism)

pupa_genome_<-filter(SRA_pupa,Organism %in% genome_list$Organism.Name) 
nrow(pupa_genome)
pupa_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

pupa_genome_list<-pupa_genome %>% select(Organism) %>% distinct(Organism)

nymph_genome<-filter(SRA_nymph,Organism %in% genome_list$Organism.Name) 
nrow(nymph_genome)
nymph_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_genome_list<-nymph_genome %>% select(Organism) %>% distinct(Organism)

larva_genome<-filter(SRA_larva,Organism %in% genome_list$Organism.Name) 
nrow(larva_genome)
larva_genome %>% select(Organism) %>% distinct(Organism) %>%nrow()

larva_genome_list<-larva_genome %>% select(Organism) %>% distinct(Organism)
```

### plot number of annotated experiments and Organisnisms

``` r
sra_samples<- c("embryo", "pupa", "nymph", "larva")

counts_org_genome<-c(embryo_genome %>% select(Organism) %>% distinct(Organism) %>%nrow(), pupa_genome %>% select(Organism) %>% distinct(Organism) %>%nrow(), nymph_genome %>% select(Organism) %>% distinct(Organism) %>% nrow(), larva_genome %>% select(Organism) %>% distinct(Organism) %>%nrow())

counts_exp_genome<-c(embryo_genome %>%nrow(), pupa_genome %>% nrow(), nymph_genome %>% nrow(), larva_genome %>%nrow())


ggplot() + geom_point(aes(x=sra_samples, y=counts_org_genome, colour  = sra_samples))  + geom_point(aes(x=sra_samples, y=counts_exp_genome, colour  = sra_samples))  + geom_point(size = 5) 
```

### filter SRA data with annotated genome which have data for embryo+ pupa and embryo+nymph/larva

``` r
embryo_genome_pupa<-filter(embryo_genome,Organism %in% pupa_genome_list$Organism) 
nrow(embryo_genome_pupa)

embryo_egg_OTHER_genome_pupa<-filter(embryo_egg_genome_hexa_OTHER,Organism %in% pupa_genome_list_hexa$Organism) 
nrow(embryo_egg_OTHER_genome_pupa)

embryo_egg_RNA_Seq_genome_pupa<-filter(embryo_egg_genome_hexa_RNA_Seq,Organism %in% pupa_genome_list_hexa$Organism) 
nrow(embryo_egg_RNA_Seq_genome_pupa)

embryo_genome_pupa %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_OTHER_genome_pupa %>% select(Organism) %>% distinct(Organism) %>%nrow()
embryo_egg_RNA_Seq_genome_pupa %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_genome_pupa_list<-embryo_genome_pupa %>% select(Organism) %>% distinct(Organism)

embryo_egg_OTHER_genome_pupa_list<-embryo_egg_OTHER_genome_pupa %>% select(Organism) %>% distinct(Organism)
embryo_egg_RNA_Seq_genome_pupa_list<-embryo_egg_RNA_Seq_genome_pupa %>% select(Organism) %>% distinct(Organism)



pupa_genome_embryo<-filter(pupa_genome,Organism %in% embryo_genome_list$Organism) 
nrow(pupa_genome_embryo)
pupa_genome_embryo %>% select(Organism) %>% distinct(Organism) %>%nrow()

pupa_genome_embryo_list<-pupa_genome_embryo %>% select(Organism) %>% distinct(Organism)

embryo_genome_pupa2 <-embryo_genome_pupa %>% mutate(stage = "embryo")  # add column named "stage" with value: embryo/pupa
pupa_genome_embryo2<- pupa_genome_embryo %>% mutate(stage = "pupa")

pupa_genome_embryo_egg_OTHER<-filter(pupa_genome_hexa,Organism %in% embryo_egg_OTHER_genome_pupa$Organism) 
nrow(pupa_genome_embryo_egg_OTHER)
pupa_genome_embryo_egg_OTHER %>% select(Organism) %>% distinct(Organism) %>%nrow()

pupa_genome_embryo_egg_OTHER_list<-pupa_genome_embryo_egg_OTHER %>% select(Organism) %>% distinct(Organism)

#
pupa_genome_embryo_egg_RNA_Seq<-filter(pupa_genome_hexa,Organism %in% embryo_egg_RNA_Seq_genome_pupa$Organism) 
nrow(pupa_genome_embryo_egg_RNA_Seq)
pupa_genome_embryo_egg_RNA_Seq %>% select(Organism) %>% distinct(Organism) %>%nrow()

pupa_genome_embryo_egg_RNA_Seq_list<-pupa_genome_embryo_egg_RNA_Seq %>% select(Organism) %>% distinct(Organism)


embryo_genome_pupa2 <-embryo_genome_pupa %>% mutate(stage = "embryo")  # add column named "stage" with value: embryo/pupa
pupa_genome_embryo2<- pupa_genome_embryo %>% mutate(stage = "pupa")

embryo_egg_OTHER_genome_pupa2 <- embryo_egg_OTHER_genome_pupa  %>% mutate(stage = "embryo")
embryo_egg_RNA_Seq_genome_pupa2 <- embryo_egg_RNA_Seq_genome_pupa  %>% mutate(stage = "embryo")

pupa_genome_embryo_egg_OTHER2 <- pupa_genome_embryo_egg_OTHER %>% mutate(stage = "pupa")
pupa_genome_embryo_egg_RNA_Seq2 <- pupa_genome_embryo_egg_RNA_Seq %>% mutate(stage = "pupa")

embryo_genome_nymph<-filter(embryo_genome,Organism %in% nymph_genome_list$Organism) 
nrow(embryo_genome_nymph)
embryo_genome_nymph %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_genome_nymph_list<-embryo_genome_nymph %>% select(Organism) %>% distinct(Organism)


embryo_egg_genome_nymph_juvenile_OTHER<-filter(embryo_egg_genome_hexa_OTHER,Organism %in% nymph_juvenile_hexa_list$Organism) 
nrow(embryo_egg_genome_nymph_juvenile_OTHER)
embryo_egg_genome_nymph_juvenile_OTHER %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_genome_nymph_juvenile_OTHER_list<-embryo_egg_genome_nymph_juvenile_OTHER %>% select(Organism) %>% distinct(Organism)

embryo_egg_genome_nymph_juvenile_RNA_Seq<-filter(embryo_egg_genome_hexa_RNA_Seq,Organism %in% nymph_juvenile_hexa_list$Organism) 
nrow(embryo_egg_genome_nymph_juvenile_RNA_Seq)
embryo_egg_genome_nymph_juvenile_RNA_Seq %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_egg_genome_nymph_juvenile_RNA_Seq_list<-embryo_egg_genome_nymph_juvenile_RNA_Seq %>% select(Organism) %>% distinct(Organism)

embryo_egg_genome_nymph_juvenile_RNA_Seq2 <-embryo_egg_genome_nymph_juvenile_RNA_Seq %>% mutate(stage = "embryo") 
embryo_egg_genome_nymph_juvenile_OTHER2 <-embryo_egg_genome_nymph_juvenile_OTHER%>% mutate(stage = "embryo") 

nymph_juvenile_genome_embryo_egg_OTHER<-filter(nymph_juvenile_hexa,Organism %in% embryo_egg_genome_hexa_OTHER_list$Organism) 
nrow(nymph_juvenile_genome_embryo_egg_OTHER)
nymph_juvenile_genome_embryo_egg_OTHER %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_juvenile_genome_embryo_egg_OTHER_list<-nymph_juvenile_genome_embryo_egg_OTHER %>% select(Organism) %>% distinct(Organism)

nymph_juvenile_genome_embryo_egg_RNA_Seq<-filter(nymph_juvenile_hexa,Organism %in% embryo_egg_genome_hexa_RNA_Seq_list$Organism) 
nrow(nymph_juvenile_genome_embryo_egg_RNA_Seq)
nymph_juvenile_genome_embryo_egg_RNA_Seq %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_juvenile_genome_embryo_egg_RNA_Seq_list<-nymph_juvenile_genome_embryo_egg_RNA_Seq %>% select(Organism) %>% distinct(Organism)

nymph_juvenile_genome_embryo_egg_RNA_Seq2<- nymph_juvenile_genome_embryo_egg_RNA_Seq %>% mutate(stage = "nympy")
nymph_juvenile_genome_embryo_egg_OTHER2<- nymph_juvenile_genome_embryo_egg_OTHER %>% mutate(stage = "nympy")

nymph_genome_embryo<-filter(nymph_genome,Organism %in% embryo_genome_list$Organism) 
nrow(nymph_genome_embryo)
nymph_genome_embryo %>% select(Organism) %>% distinct(Organism) %>%nrow()

nymph_genome_embryo_list<-nymph_genome_embryo %>% select(Organism) %>% distinct(Organism)

embryo_genome_nymph2 <-embryo_genome_nymph %>% mutate(stage = "embryo") 
nymph_genome_embryo2<- nymph_genome_embryo %>% mutate(stage = "nympy")

embryo_genome_larva<-filter(embryo_genome,Organism %in% larva_genome_list$Organism) 
nrow(embryo_genome_larva)
embryo_genome_larva %>% select(Organism) %>% distinct(Organism) %>%nrow()

embryo_genome_larva_list<-embryo_genome_larva %>% select(Organism) %>% distinct(Organism)


larva_genome_embryo<-filter(larva_genome,Organism %in% embryo_genome_list$Organism) 
nrow(larva_genome_embryo)
larva_genome_embryo %>% select(Organism) %>% distinct(Organism) %>%nrow()

larva_genome_embryo_list<-larva_genome_embryo %>% select(Organism) %>% distinct(Organism)


embryo_genome_larva2 <-embryo_genome_larva %>% mutate(stage = "embryo") 
larva_genome_embryo2<- larva_genome_embryo %>% mutate(stage = "larva")
```

### FILTERING SUMMARY

embryo_genome_pupa: runs/organism: 3106/10 pupa_genome_embryo:
runs/organism: 476/10

(before embryo_pupa_SRA - runs for embryo and pupa (13 organisms and 555
runs)

embryo_genome_nymph: runs/organism: 84/4 nymph_genome_embryo:
runs/organism: 38/4

(before embryo_nymph_SRA - runs for embryo + nymph (0 or 2 organisms and
87 runs))

embryo_genome_larva: runs/organism: 3347/18 larva_genome_embryo:
runs/organism: 1621/18

(before embryo_larva_SRA - runs for embryo + larva (14 organisms and 903
runs))

### FILTERING SUMMARY HEXA 01.2022

genome_list_hexa = 263

pupa_genome_hexa - 871 pupa_genome - 855

pupa_genome_list_hexa -45 pupa_genome_list -44

embryo_egg_genome_hexa_OTHER -5357 embryo_egg_genome_hexa_RNA_Seq -86

embryo_egg_genome_hexa_OTHER_list - 78
embryo_egg_genome_hexa_RNA_Seq_list -3

nymph_juvenile_hexa -172 nymph_juvenile_hexa_list - 17

##### 

embryo_egg_OTHER_genome_pupa - 4442 (25 organisms)
embryo_egg_RNA_Seq_genome_pupa -86

embryo_egg_OTHER_genome_pupa_list - 25 org
embryo_egg_RNA_Seq_genome_pupa_list -3 \######
pupa_genome_embryo_egg_OTHER - 653 pupa_genome_embryo_egg_OTHER_list -
25

pupa_genome_embryo_egg_RNA_Seq - 322 pupa_genome_embryo_egg_RNA_Seq- 3

embryo_egg_OTHER_genome_pupa2 embryo_egg_RNA_Seq_genome_pupa2
pupa_genome_embryo_egg_OTHER2 pupa_genome_embryo_egg_RNA_Seq2

embryo_egg_genome_nymph_juvenile_OTHER -4267
embryo_egg_genome_nymph_juvenile_OTHER_list -11

embryo_egg_genome_nymph_juvenile_RNA_Seq - 85
embryo_egg_genome_nymph_juvenile_RNA_Seq - 2

embryo_egg_genome_nymph_juvenile_RNA_Seq2
embryo_egg_genome_nymph_juvenile_OTHER2

nymph_juvenile_genome_embryo_egg_OTHER - 106
nymph_juvenile_genome_embryo_egg_OTHER_list - 11

nymph_juvenile_genome_embryo_egg_RNA_Seq - 32
nymph_juvenile_genome_embryo_egg_RNA_Seq_list -2

nymph_juvenile_genome_embryo_egg_OTHER2
nymph_juvenile_genome_embryo_egg_RNA_Seq2

##write results

``` r
library("writexl")
write_xlsx(nymph_genome_embryo2,"D:/E93_SRA_project/final/nymph_genome_embryo.xlsx")
write_xlsx(embryo_genome_nymph2,"D:/E93_SRA_project/final/embryo_genome_nymph.xlsx")

write_xlsx(embryo_genome_pupa2,"D:/E93_SRA_project/final/embryo_genome_pupa.xlsx")
write_xlsx(pupa_genome_embryo2,"D:/E93_SRA_project/final/pupa_genome_embryo.xlsx")

write_xlsx(embryo_genome_larva2,"D:/E93_SRA_project/final/embryo_genome_larva.xlsx")
write_xlsx(larva_genome_embryo2,"D:/E93_SRA_project/final/larva_genome_embryo.xlsx")

#new
write_xlsx(embryo_egg_OTHER_genome_pupa2,"HEXA_embryo_egg_OTHER_genome_pupa.xlsx")
write_xlsx(embryo_egg_RNA_Seq_genome_pupa2,"HEXA_embryo_egg_RNA_Seq_genome_pupa.xlsx")
write_xlsx(pupa_genome_embryo_egg_OTHER2,"HEXA_pupa_genome_embryo_egg_OTHER.xlsx")
write_xlsx(pupa_genome_embryo_egg_RNA_Seq2,"HEXA_pupa_genome_embryo_egg_RNA_Seq.xlsx")
write_xlsx(embryo_egg_genome_nymph_juvenile_RNA_Seq2,"HEXA_embryo_egg_genome_nymph_juvenile_RNA_Seq.xlsx")
write_xlsx(embryo_egg_genome_nymph_juvenile_OTHER2,"HEXA_embryo_egg_genome_nymph_juvenile_OTHER.xlsx")
write_xlsx(nymph_juvenile_genome_embryo_egg_OTHER2,"HEXA_nymph_juvenile_genome_embryo_egg_OTHER.xlsx")
write_xlsx(nymph_juvenile_genome_embryo_egg_RNA_Seq2,"HEXA_nymph_juvenile_genome_embryo_egg_RNA_Seq.xlsx")
```

### Plot final results

``` r
final_samples<- c("embryo_pupa",  "embryo_nymph", "embryo_larva")

counts_org<-c(embryo_genome_pupa %>% select(Organism) %>% distinct(Organism) %>%nrow(),embryo_genome_nymph %>% select(Organism) %>% distinct(Organism) %>%nrow(), embryo_genome_larva %>% select(Organism) %>% distinct(Organism) %>%nrow() )

counts_exp<-c(embryo_genome_pupa %>%nrow() + pupa_genome_embryo%>%nrow() , embryo_genome_nymph %>% nrow() + nymph_genome_embryo %>% nrow() , embryo_genome_larva %>% nrow()+larva_genome_embryo %>% nrow() )

embryo_all_nymph<-filter(SRA_embryo,Organism %in% SRA_nymph$Organism)

embryo_all_nymph_list <- embryo_all_nymph %>% select(Organism) %>% distinct(Organism)

nymph_embryo_unannotated_genome_list<- setdiff(embryo_all_nymph_list, nymph_genome_embryo_list)


write_xlsx(nymph_embryo_unannotated_genome_list,"D:/E93_SRA_project/final/nymph_embryo_unannotated_genome_list.xlsx")
```

### join data tables

``` r
embryo_genome_pupa3<-remove_empty(embryo_genome_pupa2, which = c("cols"), quiet = TRUE) # remove empty col
pupa_genome_embryo3<-remove_empty(pupa_genome_embryo2, which = c("cols"), quiet = TRUE)


pupa_final<- full_join(embryo_genome_pupa3, pupa_genome_embryo3)

nymph_genome_embryo3<-remove_empty(nymph_genome_embryo2, which = c("cols"), quiet = TRUE) # remove empty col
embryo_genome_nymph3<-remove_empty(embryo_genome_nymph2, which = c("cols"), quiet = TRUE)


nymph_final<- full_join(embryo_genome_nymph3, nymph_genome_embryo3)

embryo_genome_larva3<-remove_empty(embryo_genome_larva2, which = c("cols"), quiet = TRUE) # remove empty col
larva_genome_embryo3<-remove_empty(larva_genome_embryo2, which = c("cols"), quiet = TRUE)


larva_final<- full_join(larva_genome_embryo3, embryo_genome_larva3)

# write final results

write_xlsx(nymph_final,"D:/E93_SRA_project/final/nymph_final.xlsx")
write_xlsx(pupa_final,"D:/E93_SRA_project/final/pupa_final.xlsx")
write_xlsx(larva_final,"D:/E93_SRA_project/final/larva_final.xlsx")

# new hexa

embryo_egg_OTHER_genome_pupa3<-remove_empty(embryo_egg_OTHER_genome_pupa, which = c("cols"), quiet = TRUE)
embryo_egg_RNA_Seq_genome_pupa3<-remove_empty(embryo_egg_RNA_Seq_genome_pupa2, which = c("cols"), quiet = TRUE)

embryo_egg_genome_pupa_F<-full_join(embryo_egg_OTHER_genome_pupa3, embryo_egg_RNA_Seq_genome_pupa3)

pupa_genome_embryo_egg_OTHER3<- remove_empty(pupa_genome_embryo_egg_OTHER2, which = c("cols"), quiet = TRUE)
pupa_genome_embryo_egg_RNA_Seq3 <- remove_empty(pupa_genome_embryo_egg_RNA_Seq2, which = c("cols"), quiet = TRUE)

pupa_genome_embryo_egg_F<-full_join(pupa_genome_embryo_egg_OTHER3, pupa_genome_embryo_egg_RNA_Seq3)

PUPA_FINAL_HEXA <- full_join(embryo_egg_genome_pupa_F, pupa_genome_embryo_egg_F)

embryo_egg_genome_nymph_juvenile_RNA_Seq3 <- remove_empty(embryo_egg_genome_nymph_juvenile_RNA_Seq2, which = c("cols"), quiet = TRUE)
embryo_egg_genome_nymph_juvenile_OTHER3 <-  remove_empty(embryo_egg_genome_nymph_juvenile_OTHER2, which = c("cols"), quiet = TRUE)

embryo_egg_genome_nymph_juvenile_F <- full_join(embryo_egg_genome_nymph_juvenile_RNA_Seq3,embryo_egg_genome_nymph_juvenile_OTHER3)

nymph_juvenile_genome_embryo_egg_OTHER3 <- remove_empty(nymph_juvenile_genome_embryo_egg_OTHER2 , which = c("cols"), quiet = TRUE)
nymph_juvenile_genome_embryo_egg_RNA_Seq3  <- remove_empty(nymph_juvenile_genome_embryo_egg_RNA_Seq2 , which = c("cols"), quiet = TRUE)

nymph_juvenile_genome_embryo_egg_F<- full_join(nymph_juvenile_genome_embryo_egg_OTHER3, nymph_juvenile_genome_embryo_egg_RNA_Seq3)

NYMPH_FINAL_HEXA <- full_join(embryo_egg_genome_nymph_juvenile_F,nymph_juvenile_genome_embryo_egg_F)

holo_list<- c("Anopheles stephensi", "Apis mellifera", "Bombyx mori", "Camponotus floridanus", "Culex quinquefasciatus", "Harpegnathos saltator", "Manduca sexta", "Tribolium castaneum", "Trichomalopsis sarcophagae", "Drosophila melanogaster")


NYMPH_FINAL_HEXA2 <- NYMPH_FINAL_HEXA[ ! NYMPH_FINAL_HEXA$Organism %in% holo_list, ]

write_xlsx(NYMPH_FINAL_HEXA2,"NYMPH_FINAL_HEXA2.xlsx")
write_xlsx(PUPA_FINAL_HEXA,"PUPA_FINAL_HEXA.xlsx")
```

### Manual data validation

Metadata for each annotated species was manually reviewed to exclude
inappropriate SRA experiments, including: - non-mRNA (e.g., RIP -Seq,
EST, ncRNA-Seq) - non-embryo/pupa/nymph (juvenile) (e.g. adult, cell
culture) - treated specimens (e.g. acetone, mutants) -
non-representative body parts/tissues (e.g. wing pads, fatbody )

In some cases, where other datasets were not available for certain
species, we used data described as treated specimens (e.g., DMSO) and
non-representative body parts/tissues (e.g., fatbody ). All available
time points for embryo and pupa and the last available time point for
the nymphal/juvenile stage were used for the analysis.

### Next step –\> 2_Dowload_genomes_and_index.md
