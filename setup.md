---
title: 
---

# Setup

## Original study
This RNA-Seq lesson datasets come from a study on the model plant _Arabidopsis thaliana_ infected by a fungal pathogen called _Fusarium oxysporum_. Samples were collected from Arabidopsis root and leaf tissues from plants infected or not with the fungi after one day (1 dpi, dpi: days post-inoculation) or six days (6 dpi). 

The study was conducted by [Lyons et al. in 2015 and was published in PLoS One](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4388846).  

## Files
We will use two data files for this lesson:
- **Counts:** A `counts.tsv` dataframe of the sample raw counts. `tsv` stands for Tab Separated Values and displays data in tabulated separated columns.
- **Experimental design**: the `design.tsv` dataframe indicates the correspondence between samples and experimental conditions (e.g. control, treated).   

These two files were obtained by running the `v0.1.1` version of a [RNA-Seq bioinformatic pipeline](https://github.com/KoesGroup/Snakemake_hisat-DESeq/blob/master/README.md) on the mRNA-Seq sequencing files from Lyons et al. (2015): https://www.ncbi.nlm.nih.gov/bioproject/272563.


## Softwares

> ## What you need to install.
>
> Before the training, please make sure you have done the following: 
>
> 1. Download and install **up-to-date versions** of:
>    - R: [https://cloud.r-project.org](https://cloud.r-project.org).
>    - RStudio: [http://www.rstudio.com/download](http://www.rstudio.com/download). 
>    - The DESeq2 package: [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
> 2. Read the workshop [Code of Conduct](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html) to make sure this workshop stays welcoming for everybody.
> 3. Get comfortable: if you're not in a physical workshop, be set up with two screens if possible. You will be following along in RStudio on your own computer while also following this tutorial on your own.
> More instructions are available on the workshop website in the **Setup** section.
{: .prereq}

