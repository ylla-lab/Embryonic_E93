## Load libraries

``` r
library(DESeq2)
library(dplyr)
library(ggplot2)
library(Rsubread)
library(ggrepel)
library(edgeR)
library(limma)
library(tidyr)
library(readr)
library(sva)
library(ggthemes)
library(readxl)
```

# RNA-Seq analysis of E93-depleted and control embryos from Blatella germanica

## Get Table of counts

-   With FeatureCounts

``` r
BAMfiles1<- list.files(path="../../sam_bam/clean",pattern = "*.bam$",full.names=TRUE)
BAMfiles2<- list.files(path="../../sam_bam/clean_new",pattern = "*.bam$",full.names=TRUE)
BAMfiles<-append(BAMfiles1,BAMfiles2)
```

-   Estimate counts with FeatureCounts

``` r
annotation<-"../../genome/BGER_v0.6.2.3.1.gff"

dsE03exp_counts<-featureCounts(files=BAMfiles, isPairedEnd=TRUE, 
                              requireBothEndsMapped=FALSE,
                              allowMultiOverlap=FALSE, 
                              countMultiMappingReads=TRUE, 
                              annot.ext=annotation,
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType="gene",
                              GTF.attrType="ID",
                              useMetaFeatures=TRUE,
                              nthreads=10)
```

## Table of counts

``` r
dsE03exp_counts_df<-dsE03exp_counts$counts
colnames(dsE03exp_counts_df)<-paste0("Sample_",sapply(strsplit(colnames(dsE03exp_counts_df),"_clean"), `[`, 1))

write.csv(dsE03exp_counts_df,"E93_RNASeq_counts.csv", row.names = T)
```

``` r
Meatdata_ds3exp<- read_delim("samples_metadata.txt", 
     delim = "\t", escape_double = FALSE, 
     trim_ws = TRUE)

Meatdata_ds3exp
```

-   As DF for DESeq2

``` r
Meatdata_ds3exp_DF<-as.data.frame(Meatdata_ds3exp)
Meatdata_ds3exp_DF$SampleName<-paste0("Sample_",Meatdata_ds3exp_DF$sample)
rownames(Meatdata_ds3exp_DF)<-Meatdata_ds3exp_DF$SampleName
Meatdata_ds3exp_DF
```

## DESeq2

``` r
all(rownames(Meatdata_ds3exp_DF) == colnames(dsE03exp_counts_df))

Meatdata_ds3exp_DF_nz<-dsE03exp_counts_df[which(rowSums(dsE03exp_counts_df)>0),]

dds_dsE93 <- DESeq2::DESeqDataSetFromMatrix(countData = dsE03exp_counts_df,
                              colData = Meatdata_ds3exp_DF,
                              design= ~condition)
```

### VST transformation raw counts

``` r
dsE93exp_VST<-assay(varianceStabilizingTransformation(dds_dsE93, blind = TRUE)) 
```

-   VST long format and metadata

``` r
dsE93exp_VST_long<-dsE93exp_VST %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "GeneID") %>% 
  pivot_longer(cols=colnames(dsE93exp_VST),names_to="Sample", values_to="VST") %>% 
   dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 
```

### PCA

``` r
dsE93exp_VST_toPCA_t<-t(dsE93exp_VST)

dsE93exp_VST_toPCA<-as.data.frame(dsE93exp_VST_toPCA_t) %>% 
    tibble::rownames_to_column(var ="Sample") %>% 
    dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 

rownames(dsE93exp_VST_toPCA)<-dsE93exp_VST_toPCA$Sample


pca_VST_dsE93exp<-prcomp(dsE93exp_VST_toPCA[,2:nrow(dsE93exp_VST)+1] )



var_explained <- pca_VST_dsE93exp$sdev^2/sum(pca_VST_dsE93exp$sdev^2)
var_explained[1:5]

pca_VST_dsE93exp_plot<-pca_VST_dsE93exp$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       title = "PCA raw",color ="Sample") +
 geom_text_repel(label = dsE93exp_VST_toPCA$Sample, size=4) +
 geom_point(aes(color=dsE93exp_VST_toPCA$condition), size = 4)+
 theme_bw(base_size = 12)+
 scale_color_tableau()


pca_VST_dsE93exp_plot

svg("pca_VST_dsE93exp_plot.svg")
plot(pca_VST_dsE93exp_plot)

dev.off()
```

## Calculate SV and remove batch effect on counts (for PCA)

-   calculate SV

``` r
dds <- estimateSizeFactors(dds_dsE93)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

n.sv = num.sv(dsE03exp_counts_df,mod,method="leek") #  this comand give 10, but we use 3


svseq <- svaseq(dat, mod, mod0, n.sv=3)

ddssva <-dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

cov = cbind(svseq$sv[,1:3])
```

-   remove batch on vst

``` r
dsE93exp_VST_batch_corrected <- removeBatchEffect(dsE93exp_VST, covariates = cov) 
```

-   VST long format and metadata

``` r
dsE93exp_VST_batch_corrected_long<-dsE93exp_VST_batch_corrected  %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "GeneID") %>% 
  pivot_longer(cols=colnames(dsE93exp_VST_batch_corrected ),names_to="Sample", values_to="VST") %>% 
   dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 
```

### PCA

``` r
dsE93exp_VST_batch_corrected_toPCA_t<-t(dsE93exp_VST_batch_corrected)

dsE93exp_VST_batch_corrected_toPCA<-as.data.frame(dsE93exp_VST_batch_corrected_toPCA_t) %>% 
    tibble::rownames_to_column(var ="Sample") %>% 
    dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 

rownames(dsE93exp_VST_batch_corrected_toPCA)<-dsE93exp_VST_batch_corrected_toPCA$Sample


pca_VST_dsE93exp_batch_corrected<-prcomp(dsE93exp_VST_batch_corrected_toPCA[,2:nrow(dsE93exp_VST_batch_corrected)+1] )

var_explained2 <- pca_VST_dsE93exp_batch_corrected$sdev^2/sum(pca_VST_dsE93exp_batch_corrected$sdev^2)
var_explained2[1:5]
```

``` r
pca_VST_dsE93exp_batch_corrected_plot_a<-pca_VST_dsE93exp_batch_corrected$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  labs(x=paste0("PC1: ",round(var_explained2[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained2[2]*100,1),"%"),
       title = "PCA afrer batch correction",color ="Sample") +
 geom_text_repel(label = dsE93exp_VST_batch_corrected_toPCA$Sample, size=4) +
 geom_point(aes(color=dsE93exp_VST_batch_corrected_toPCA$condition), size = 4)+
 theme_bw(base_size = 12)+
 scale_color_tableau()


pca_VST_dsE93exp_batch_corrected_plot_a

#svg("pca_VST_dsE93exp_batch_corrected_plot_all.svg")
#plot(pca_VST_dsE93exp_batch_corrected_plot_a)

#dev.off()
```

``` r
library(ggpubr)
    
plots_supl_PCA<- ggarrange(pca_VST_dsE93exp_plot, pca_VST_dsE93exp_batch_corrected_plot_a,nrow = 1, ncol=2,labels = c("A", "B") , common.legend = T, legend = "bottom" )
plots_supl_PCA



svg("plots_supl_PCA2.svg")
plot(plots_supl_PCA)
dev.off()
```

## Remove sample T7 based on PCA

``` r
dsE03exp_counts_df<-dplyr::select(as.data.frame(dsE03exp_counts_df), Sample_1, Sample_2, Sample_3, Sample_6,Sample_8, Sample_9, Sample_C9, Sample_C11, Sample_C12, Sample_T6, Sample_T14)
```

## Load Metdata

-   Meatdata

``` r
Meatdata_ds3exp <- read_delim("samples_metadata_NO_T7.txt", 
     delim = "\t", escape_double = FALSE, 
     trim_ws = TRUE)

Meatdata_ds3exp
```

-   As DF for DESeq2

``` r
Meatdata_ds3exp_DF<-as.data.frame(Meatdata_ds3exp)
Meatdata_ds3exp_DF$SampleName<-paste0("Sample_",Meatdata_ds3exp_DF$sample)
rownames(Meatdata_ds3exp_DF)<-Meatdata_ds3exp_DF$SampleName
Meatdata_ds3exp_DF
```

## DESeq2

``` r
all(rownames(Meatdata_ds3exp_DF) == colnames(dsE03exp_counts_df))

Meatdata_ds3exp_DF_nz<-dsE03exp_counts_df[which(rowSums(dsE03exp_counts_df)>0),]

dds_dsE93 <- DESeq2::DESeqDataSetFromMatrix(countData = dsE03exp_counts_df,
                              colData = Meatdata_ds3exp_DF,
                              design= ~condition)
```

## Calculate SV and remove batch effect on counts (for PCA)

-   calculate SV

``` r
dds <- estimateSizeFactors(dds_dsE93)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

n.sv = num.sv(dsE03exp_counts_df,mod,method="leek") #  this comand give 9, but we use 3


svseq <- svaseq(dat, mod, mod0, n.sv=3)

ddssva <-dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

cov = cbind(svseq$sv[,1:3])
```

-   VST transformation

``` r
dsE93exp_VST<-assay(varianceStabilizingTransformation(dds_dsE93, blind = TRUE)) 
```

-   remove batch on vst

``` r
dsE93exp_VST_batch_corrected <- removeBatchEffect(dsE93exp_VST, covariates = cov) 
```

-   VST long format and metadata

``` r
dsE93exp_VST_batch_corrected_long<-dsE93exp_VST_batch_corrected  %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "GeneID") %>% 
  pivot_longer(cols=colnames(dsE93exp_VST_batch_corrected ),names_to="Sample", values_to="VST") %>% 
   dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 
```

### PCA

``` r
dsE93exp_VST_batch_corrected_toPCA_t<-t(dsE93exp_VST_batch_corrected)

dsE93exp_VST_batch_corrected_toPCA<-as.data.frame(dsE93exp_VST_batch_corrected_toPCA_t) %>% 
    tibble::rownames_to_column(var ="Sample") %>% 
    dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 

rownames(dsE93exp_VST_batch_corrected_toPCA)<-dsE93exp_VST_batch_corrected_toPCA$Sample


pca_VST_dsE93exp_batch_corrected<-prcomp(dsE93exp_VST_batch_corrected_toPCA[,2:nrow(dsE93exp_VST_batch_corrected)+1] )

var_explained2 <- pca_VST_dsE93exp_batch_corrected$sdev^2/sum(pca_VST_dsE93exp_batch_corrected$sdev^2)
var_explained2[1:5]
```

``` r
pca_VST_dsE93exp_batch_corrected_plot_no_T7<-pca_VST_dsE93exp_batch_corrected$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + 
  labs(x=paste0("PC1: ",round(var_explained2[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained2[2]*100,1),"%"),
       title = "PCA afrer batch correction",color ="Sample") +
 geom_text_repel(label = dsE93exp_VST_batch_corrected_toPCA$Sample, size=4) +
 geom_point(aes(color=dsE93exp_VST_batch_corrected_toPCA$condition), size = 4)+
 theme_bw(base_size = 12)+
 scale_color_tableau()


pca_VST_dsE93exp_batch_corrected_plot_no_T7

svg("pca_VST_dsE93exp_batch_corrected_plot_no_T7.svg")
plot(pca_VST_dsE93exp_batch_corrected_plot_no_T7)

dev.off()
```

## Remove batch for DEA

``` r
dds <- estimateSizeFactors(dds_dsE93)
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

n.sv = num.sv(dsE03exp_counts_df,mod,method="leek") #  this comand give 9, but we use 3

svseq2 <- svaseq(dat, mod, mod0, n.sv=3)

ddssva2 <-dds
ddssva2$SV1 <- svseq2$sv[,1]
ddssva2$SV2 <- svseq2$sv[,2]
ddssva2$SV3 <- svseq2$sv[,3]
```

## DEA

``` r
design(ddssva2) <- ~ SV1 + SV2 + SV3 + condition

ddssva2 <- DESeq(ddssva2)

resultsNames(ddssva2)

ddssva2_result<- results(ddssva2, contrast = c("condition", "dsE93","control"))
```

``` r
ddssva2_result_tbl<-ddssva2_result %>%as.data.frame() %>%   
  tibble::rownames_to_column(var="GeneID") %>%as_tibble() %>%  
  dplyr::arrange(padj) %>% 
  dplyr::mutate(Effect = dplyr::if_else( log2FoldChange>0, "dsE93_Up","dsE93_Down" ))
```

### VST

``` r
ddssva2_VST<-assay(varianceStabilizingTransformation(ddssva2, blind = FALSE))
```

-   VST long format and metadata

``` r
ddssva2_VST_long<-ddssva2_VST %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "GeneID") %>% 
  pivot_longer(cols=colnames(ddssva2_VST),names_to="Sample", values_to="VST") %>% 
   dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 
```

### Add orthology info

``` r
Orthology_df<-read.csv("Orthologs_Names.csv")
```

``` r
ddssva2_result_tbl_annot<-ddssva2_result_tbl %>% left_join(Orthology_df, by=c("GeneID"="Bger"))
```

-   summarize

``` r
ddssva2_result_tbl %>% 
  filter(padj<0.05) %>% 
  group_by(Effect) %>% 
  dplyr::summarise(n())
```

``` r
MA_plot<-DESeq2 :: plotMA(ddssva2_result , colNonSig = "grey")

svg("MA_plot.svg")
plot(MA_plot)

dev.off()
```

-   write results tables

``` r
DEG_0.05<-ddssva2_result_tbl_annot %>% 
  filter(padj<0.05) %>% arrange(desc(abs(log2FoldChange)))

write.csv(DEG_0.05, file = "DEG_0.05.csv", row.names = T)
write.csv(ddssva2_result_tbl_annot, file = "DEA_results.csv", row.names = T)
```

# Draw expression plots on VST bath corrected

-   VST long format and metadata

``` r
dsE93exp_VST_batch_corrected_long<-dsE93exp_VST_batch_corrected %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "GeneID") %>% 
  pivot_longer(cols=colnames(dsE93exp_VST_batch_corrected),names_to="Sample", values_to="VST") %>% 
   dplyr::left_join(Meatdata_ds3exp_DF, by=c("Sample"="SampleName")) 
```

-   draw expression plots on VST

``` r
dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="E93") %>% 
  ggplot(., aes(x=Sample  ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "E93 expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_01853") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Kruppel expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()
    

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Scaffold533:124865-197312") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Orthodenticle/Ocelliless expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Caudal") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Caudal expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()



dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_01219") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "ActinC5 expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_18070") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Even skipped expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_23144") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Nanos expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()


dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Zelda") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Zelda expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_03424") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Knirps expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()


dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_14305") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Hunchback expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_19904") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Tailless - NO expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_26281") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Hairy expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()


dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_19904") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Fushi tarazu - NO expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_05999") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Runt  expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_14451") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Odd-paired  expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()


dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_03163") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Sloppy paired 2 expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()

dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_00720") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Engrailed expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()


dsE93exp_VST_batch_corrected_long%>%
  filter(GeneID=="Bger_05081") %>% 
  ggplot(., aes(x=Sample   ,color=condition, y=VST)) +
    geom_point( position = position_jitterdodge(), size = 8)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 35, hjust = .8, size=10))+
    labs(title = "Gooseberry-neuro expression",y="Normalized counts (corrected vst)",fill ="Legend", shape="Shape")+ #+guides(color=FALSE)
    scale_color_tableau()
```

# ClusterProfiler GSEA

-   filtered with padj\<1 ( to filter out NA)

``` r
filtered_padj<-ddssva2_result_tbl_annot %>% 
  filter(padj<1) %>% 
  arrange(desc(abs(log2FoldChange)))
```

-   read data

``` r
# log2 fold change 
original_gene_list <- filtered_padj$log2FoldChange

# name the vector
names(original_gene_list) <- filtered_padj$Dmel

# omit NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order
gene_list = sort(gene_list, decreasing = TRUE)

ego2BP <- gseGO(gene     = gene_list,
               OrgDb    = org.Dm.eg.db,
              keyType       = 'ENSEMBL',
               ont      = "BP",  # biological process
               pAdjustMethod = "BH",
              nPermSimple = 10000,
               pvalueCutoff  = 0.05)

ego2BP_df<-as.data.frame(ego2BP)

ego2MF <- gseGO(gene     = gene_list,
               OrgDb    = org.Dm.eg.db,
              keyType       = 'ENSEMBL',
               ont      = "MF",  # molecular function
               pAdjustMethod = "BH",
              nPermSimple = 10000,
               pvalueCutoff  = 0.05)


ego2MF_df<-as.data.frame(ego2MF)

ego2CC <- gseGO(gene     = gene_list,
               OrgDb    = org.Dm.eg.db,
              keyType       = 'ENSEMBL',
               ont      = "CC",  # molecular function
               pAdjustMethod = "BH",
              nPermSimple = 10000,
               pvalueCutoff  = 0.05)


ego2CC_df<-as.data.frame(ego2CC)
```

-   write results table BP and MF

``` r
write.csv(ego2BP_df, file = "BP_padj_filtered.csv", row.names = T)
write.csv(ego2MF_df, file = "MF_padj_filtered.csv", row.names = T)
write.csv(ego2CC_df, file = "CC_padj_filtered.csv", row.names = T)
```

-   Dotplot

``` r
require(DOSE)
dotplot(ego2BP, showCategory=15, split=".sign",font.size =7, title="Enrichment top 15 BP") + facet_grid(.~.sign)

dotplot(ego2MF, showCategory=100, split=".sign", font.size =7, title="Enrichment MF") + facet_grid(.~.sign)

dotplot(ego2CC, showCategory=15, split=".sign", font.size =7, title="Enrichment CC") + facet_grid(.~.sign)


plot_BP_15<-dotplot(ego2BP, showCategory=20, split=".sign",font.size =7.5, title="Enrichment top 20 BP") + facet_grid(.~.sign)

svg("plot_BP_15.svg")
plot(plot_BP_15)

dev.off()

plot_MF<-dotplot(ego2MF, showCategory=100, split=".sign", font.size =7.5, title="Enrichment MF") + facet_grid(.~.sign)

svg("plot_MF.svg")
plot(plot_MF)

dev.off()

plot_CC<-dotplot(ego2MF, showCategory=15, split=".sign", font.size =7.5, title="Enrichment top 20 CC ") + facet_grid(.~.sign)

svg("plot_CC_15.svg")
plot(plot_CC)

dev.off()

library(ggpubr)
    
plots_supl_GSEA1<- ggarrange(plot_BP_15,nrow = 1, ncol=1,labels = c("A") , common.legend = F, legend = "right" )
plots_supl_GSEA1
plots_supl_GSEA2<- ggarrange(plot_CC,nrow = 1, ncol=1,labels = c("B") , common.legend = F, legend = "right" )
plots_supl_GSEA2
plots_supl_GSEA3<- ggarrange(plot_MF,nrow = 1, ncol=1,labels = c("C") , common.legend = F, legend = "right" )
plots_supl_GSEA3

svg("plots_supl_GSEA2.svg")
plot(plots_supl_GSEA2)
dev.off()
```

# KEGG

``` r
# Convert gene IDs for gseKEGG function

ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Dm.eg.db)

 # remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 with mapped genes
df2 =  filtered_padj[ filtered_padj$Dmel %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
```

``` r
kegg_organism = "dme"
kk2_sel <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kk2_sel_df<-as.data.frame(kk2_sel)
```

``` r
write.csv(kk2_sel_df, file = "KEGG.csv", row.names = T)
```

``` r
dotplot(kk2_sel, showCategory = 100, title = "Enriched Pathways" , split=".sign",  font.size =6) + facet_grid(.~.sign)

plot_KEGG<-dotplot(kk2_sel, showCategory = 100, title = "Enriched Pathways" , split=".sign",  font.size =7.5) + facet_grid(.~.sign)

svg("plot_KEGG.svg")
plot(plot_KEGG)

dev.off()
```

``` r
plot1_dsE93<- ggarrange(pca_VST_dsE93exp_batch_corrected_plot_no_T7, MA_plot ,nrow = 2, ncol=1,labels = c("A","B"), common.legend = F, legend = "right" )
plot1_dsE93


plot1_dsE932<- ggarrange(plot1_dsE93, plot_KEGG ,nrow = 1, ncol=2,labels = c("","C"), common.legend = F, legend = "right" )
plot1_dsE932


svg("plot1_dsE932.svg")
plot(plot1_dsE932)
dev.off()
```

## Session Info

``` r
sessionInfo()
```
