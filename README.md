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


## Dataset
The data used for analysis is from the study, “EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival” (Fu et al. 2015).

This study examines the expression profiles of basal stem-cell enriched cells  and committed luminal cells in the mammary gland of virgin, pregnant and lactating mice.

GEO Accession ID - [GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450)




## Results
The differentially expressed genes were obtained. 
Various plots like heatmap, mean-variance plot, MA plot, and volcano plot were made to analyse the data better.


### Library Size of the Samples
[![Picture4.png](https://i.postimg.cc/DykqwfhK/Picture4.png)](https://postimg.cc/0rn6tsrt)

### HeatMap
[![Picture3.png](https://i.postimg.cc/v8zQc3FH/Picture3.png)](https://postimg.cc/w31Kr5jS)

### Mean-Variance Plot
[![Picture1.png](https://i.postimg.cc/50QCfTk0/Picture1.png)](https://postimg.cc/ygBdPQK4)

### MA Plot & Volcano Plot
[![Picture2.png](https://i.postimg.cc/LXFHZW8h/Picture2.png)](https://postimg.cc/7CXvFm7y)





The dataset used here is present in `Datasets` folder.

