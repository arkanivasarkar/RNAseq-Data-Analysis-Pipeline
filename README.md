# RNAseq Data Analysis Pipeline

This repository contains the workflow of RNAseq data analysis.

This pipeline performs the following tasks:

- Reading data
- align reads of each sample in a run against reference genome
= perform quality control on generated BAM files 
- count reads in features 
- normalize read counts
- Filtering lowly expressed genes
- perform DE analysis 


### Dataset
The data used for analysis is from the study, “EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival” (Fu et al. 2015).

This study examines the expression profiles of basal stem-cell enriched cells  and committed luminal cells in the mammary gland of virgin, pregnant and lactating mice.

GEO Accession ID - [GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)




## Results
The differentially expressed genes were obtained. 
Various plots like heatmap, mean-variance plot, MA plot, and volcano plot were made to analyse the data better.


# Library Size of the Samples
![alt text](https://postimg.cc/0rn6tsrt) 

# HeatMap
![alt text](https://i.ibb.co/BKr7cVF/Picture1.png) 

# Mean-Variance Plot
![alt text](https://postimg.cc/ygBdPQK4) 

# MA Plot & Volcano Plot
![alt text](https://i.ibb.co/BKr7cVF/Picture1.png) 





The dataset used here is present in `Datasets` folder.

