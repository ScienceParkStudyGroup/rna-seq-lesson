---
layout: lesson
authors: Marc Galland, Tijs Bliek, Ken Kraaijeveld
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---
# Welcome!

## Overview 

In this training you will learn the basics of a typical RNA-Sequencing experiment. It is going to be fun and empowering! You will discover how total RNA are converted to short sequences called "reads" that can in turn be used to get insights into gene expression. Through careful experimental design, these gene expression information can yield new research avenues and answer crucial questions. 

We will practice learning three main things all at the same time: 
1. The purpose and good practices related to an RNA-Seq experiment. 
2. Making plots to describe the output of a RNA-Seq workflow.
3. How to perform a differential expression analysis and create some insightful plots for one comparison.

We will use mostly R and its companion RStudio to perform our RNA-Seq analyses and visualisations. 

Depending on the level of participants, the bioinformatic part might be performed (QC of fastq files, genome alignment, counting, etc.)  

This lesson will introduce you to the basics of gene expression analysis using RNA-Seq (short for RNA sequencing). Due to the considerable progress and constant decreasing costs of RNA-Seq, this technique has became a standard

## Main learning objectives
- Identify good practices when designing a RNA-Seq experiment.
- Memorize the steps of a complete RNA-Seq experiment: from sequencing to analysis.
- Perform a QC of your experiment through Principal Component Analysis (PCA) and sample clustering. 
- Execute a differential gene expression analysis using R and the `DESeq2` package.  
- Be able to create the key plots: volcano plot, heatmap and clustering of differentially expressed genes.  


Before you begin, be sure you are all set up: see the prerequisites (see below) and in the setup section. 
> ## Before you start
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


{% include links.md %}


----



