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
- "Tidy data principles are essential to increase data analysis efficiency and code readability."
- "Using R and RStudio, it becomes easier to implement good practices in data analysis."
- "I can make my workflow more reproducible and collaborative by using git and Github."
---

## Table of Contents
1. [Overview](#overview)
2. [What you will learn](#what-you-will-learn)
3. [Dataset used](#dataset-used)


## Overview 

Welcome!

In this training you will learn the basics of a typical RNA-Sequencing experiment. It is going to be fun and empowering! You will discover how total RNA are converted to short sequences called "reads" that can in turn be used to get insights into gene expression. Through careful experimental design, these gene expression information can yield new research avenues and answer crucial questions. 

We will practice learning three main things all at the same time: 
1. The purpose and good practices related to an RNA-Seq experiment. 
2. Making plots to describe the output of a RNA-Seq workflow.
3. How to perform a differential expression analysis and create some insightful plots for one comparison.

We will use mostly R and its companion RStudio to perform our RNA-Seq analyses and visualisations. 

Depending on the level of participants, the bioinformatic part might be performed (QC of fastq files, genome alignment, counting, etc.)  

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

Arabidopsis plants 


![](../img/experimental_design.png)


The original sequencing files can be found on the [Array Express database of the European Bioinformatic Institute](https://www.ebi.ac.uk/arrayexpress) by searching for the dataset accession number __E‐MTAB‐4683__.

## Credits

### Dataset
The original RNA-Seq dataset used comes from Vogel et al. 2016:  https://doi.org/10.1111/nph.14036.  

### Teaching materials
This lesson has been formatted according to the [Carpentries Foundation](https://carpentries.org/) lesson template and following their recommendations on how to teach researchers good practices in programming and data analysis.   

This material builds from a lot of fantastic materials developed by others in the open data science community. Most of the content derives from the [Harvard Chan Bioinformatics Core](https://github.com/hbctraining) which are greatly acknowledge for the quality of their teaching materials.

It also pulls from the following resources, which are highly recommended for further learning and as resources later on. Specific lessons will also cite more resources.

- [R for Data Science](http://r4ds.had.co.nz/) by Hadley Wickham and Garrett Grolemund
- [STAT 545](http://stat545.com/) by Jenny Bryan
- [Happy Git with R](http://happygitwithr.com) by Jenny Bryan
- [Software Carpentry](https://software-carpentry.org/lessons/) by the Carpentries

{% include links.md %}

<!---
## Motivation 


More often than not, there are more than one way to do things. I'm going to focus mostly on what I have ended up using day-to-day; I try to incorporate better practices as I come upon them but that's not always the case. RStudio has some built-in redundancy too that I'll try to show you so that you can approach things in different ways and ease in.

- based on literature: best and good enough practices
- also based on our team's experience of how to do better science in less time




## Collaboration

Everything we learn today is to going to help you collaborate with your most important collaborator — YOU. Science is collaborative, starting with Future You, your current collaborators, and anyone wanting to build off your science later on. 

## Reproducibility

- record of your analyses. 
- rerun them!
- modify them, maybe change a threshold, try a different coefficient, etc, maybe today
- modify them, make a new figure, in 6 months! 

## Mindset

New but will become increasingly familiar. We’ll start you off with some momentum, like if you were going to learn to ride a bike or ...

Expect that there is a way to do what you want to do

- stop confounding data science with your science. Expect that someone has had your problem before or done what you want to do. 


If you plan to program mostly in one particular language on a single platform (such as Mac or Windows), you might try an integrated development environment (IDE). IDEs integrate text editing, syntax highlighting, version control, help, build tools, and debugging in one interface, simplifying development. 

http://r-bio.github.io/intro-git-rstudio/

## Data science is a discipline

It has theories, methods, and tools. 

Tidyverse and Hadley’s graphic. Tidy data.

Going to teach you how to think differently, get into some of the theory but in the context of hands-on work.


--->
