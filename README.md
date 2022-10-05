# Reduction of embryonic E93 expression as a key factor for the evolution of insect metamorphosis

--------------

![pubstatus](https://img.shields.io/badge/Article_Status:-Sent_to_journal-orange)

--------------


- Article *"[Reduction of embryonic E93 expression as a key factor for the evolution of insect metamorphosis](https://www.biorxiv.org/content/10.1101/2022.10.04.510826v1)"*
- Link to preprint in bioRxiv: https://www.biorxiv.org/content/10.1101/2022.10.04.510826v1
- Scripts author: Gabriela Machaj, 2022
- [Laboratory of Bioinformatics and Genome Biology](https://ylla-lab.github.io/), Jagiellonian University

--------------

This repository contains the scripts written by Dr. Machaj for analyzing the data and obtaining the results publsihed at *"Reduction of embryonic E93 expression as a key factor for the evolution of insect metamorphosis"*. 

Scripts are organized in two main analysis:

1. [RNA-seq analysis of *B. germanica* libraries under E93 knockdown compared to control](https://github.com/ylla-lab/Embryonic_E93/tree/master/RNAseq_Bger_dsE93)
	1. [Prepare data, QC, and mapping](RNAseq_Bger_dsE93/1_Data_download_to_mapping.md )
	2. [Table of counts, DEA, GSEA, etc.](https://htmlpreview.github.io/?https://raw.githubusercontent.com/ylla-lab/Embryonic_E93/master/RNAseq_Bger_dsE93/2_E93_R_read_count_to_GSEA.html)
2. [Comparison of E93 expression in embryo vs pre-adult across hexapods](https://github.com/ylla-lab/Embryonic_E93/tree/master/RNAseq_multipleSpp_E93exp)
	1. [Identify appropriate datasets in NCBI-SRA](https://github.com/ylla-lab/Embryonic_E93/blob/master/RNAseq_multipleSpp_E93exp/1_Filtering-SRA-data.md)
	2. [Download genomes form NCBI](https://github.com/ylla-lab/Embryonic_E93/blob/master/RNAseq_multipleSpp_E93exp/2_Download_genomes_and_index.md)
	3. [Download and map RNA-seq data from NCBI-SRA](https://github.com/ylla-lab/Embryonic_E93/blob/master/RNAseq_multipleSpp_E93exp/3_Download_SRA_data_filter_gality_control_mapping.md)
	4. [Generate tables of counts](https://github.com/ylla-lab/Embryonic_E93/blob/master/RNAseq_multipleSpp_E93exp/4_Generate_table_of_counts.Rmd)
	5. [Calculate E93 expression ratio in pre-adult vs embryo](https://htmlpreview.github.io/?https://raw.githubusercontent.com/ylla-lab/Embryonic_E93/master/RNAseq_multipleSpp_E93exp/5_SRA_insects_E93_ratios.html)



