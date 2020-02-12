---
title: "Introduction"
teaching: 30
exercises: 0
questions:
- "What can I learn by doing this RNA-Seq lesson?"
- "What are the tools that I will be using?"
- "What are the tidy data principles?"
- "What is working in a more open way beneficial?"
objectives:
- "Discover in what consists a typical RNA-Seq experiment."
- "Know what results a RNA-Seq experiment can (and cannot) yield."
keypoints:
- "An RNA-Seq experiment is also a normal experiment (control, treated, replication, etc.)."
- A canonical RNA-Seq experiment consists in RNA library preparation followed by bioinformatic analyses."
- "RNA-Seq yields a snapshot of individual gene expression levels (count table)."
- "Upon completion of the bioinformatic steps, the analysis of RNA-Seq results can be done using the DESeq2 R package."
---

## Table of Contents
1. [Overview](#overview)
2. [What you will learn](#what-you-will-learn)
3. [Dataset used](#dataset-used)


## What you will learn

1. **What are the important things to know before doing an RNA-Seq experiment** 
    - When should you perform a RNA-Seq experiment?  
    - RNA-Seq experiments have to comply with good experimental design practices just like any experiment.
    - What are biological replicates and why are they important?
2. **How can I assess the quality of my RNA-Seq sequencing results?**
    - FastQC.
    - PCA plot.
    - Sample clustering.
3. **How do I perform a differential expression analysis on RNA-Seq results using R?**
    - Raw and scaled counts: why do you need to scale counts between samples?
    - What are the gene expression units I need to know: RPKM, FPKM, TPM.
    - What are robust scaling/normalisation methods?
    - How does the DESeq method works?
4. **What are the plots that I can make from the differential analysis results?**
    - Heatmap coupled with gene and sample clustering.
    - Volcano plot.

## Dataset used 

We will make use of a published experimental dataset from a study made on the small model plant _Arabidopsis thaliana_ by [Vogel et al. (2016)](https://doi.org/10.1111/nph.14036). This study compares the response of 4 weeks old plantlets to different bacteria that live on the leaves of different plant species:
- A known foliar pathogen called _Pseudomonas syringae_ strain DC3000. 
- A commensal ("neutral") bacteria called _Methylobacterium extorquens_ strain PA1.
- A commensal ("neutral") bacteria called _Sphingomonas melonis_ strain Fr1. 

<img src="../img/experimental_design.png" width="600px" alt="experimental design" >


The original sequencing files can be found on the [Array Express database of the European Bioinformatic Institute](https://www.ebi.ac.uk/arrayexpress) by searching for the dataset accession number __E‐MTAB‐4683__.

## Credits

### Dataset
The original RNA-Seq dataset used comes from Vogel et al. 2016:  https://doi.org/10.1111/nph.14036.  

### Teaching materials
This lesson has been formatted according to the [Carpentries Foundation](https://carpentries.org/) lesson template and following their recommendations on how to teach researchers good practices in programming and data analysis.   

This material builds from a lot of fantastic materials developed by others in the open data science community. Most of the content derives from the [Harvard Chan Bioinformatics Core](https://github.com/hbctraining) which are greatly acknowledge for the quality of their teaching materials.

{% include links.md %}
