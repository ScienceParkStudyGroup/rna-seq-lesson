---
layout: lesson
authors: Marc Galland, Tijs Bliek, Ken Kraaijeveld
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---

## Welcome!

In this training you will learn the basics of a typical RNA-Sequencing experiment. It is going to be fun and empowering! You will discover how total RNA are converted to short sequences called "reads" that can in turn be used to get insights into gene expression. Through careful experimental design, these gene expression information can yield new research avenues and answer crucial questions. 

We will use mostly R and its companion RStudio to perform our RNA-Seq analyses and visualisations. 

Depending on the level of participants, the bioinformatic part might be performed (QC of fastq files, genome alignment, counting, etc.)

Before you begin, be sure you are all set up (see below). For complete information, see the Setup section. 

We will use mostly R and its companion RStudio to perform our RNA-Seq analyses and visualisations. 

Depending on the level of participants, the bioinformatic part might be performed (QC of fastq files, genome alignment, counting, etc.)  

This lesson will introduce you to the basics of gene expression analysis using RNA-Seq (short for RNA sequencing). Due to the considerable progress and constant decreasing costs of RNA-Seq, this technique has became a standard

## Main learning objectives
- Identify good practices when designing a RNA-Seq experiment.
- Memorize the steps of a complete RNA-Seq experiment: from sequencing to analysis.
- Perform a QC of your experiment through Principal Component Analysis (PCA) and sample clustering. 
- Execute a differential gene expression analysis using R and the `DESeq2` package.  
- Be able to create key plots: volcano plot, heatmap and clustering of differentially expressed genes.  
- Provide a biological interpretation to differentially expressed genes through ORA/GSEA analyses and data integration.

> ## Before you start
>
> Before the training, please make sure you have done the following: 
>
> 1. Consult what you need to do in the lesson [Setup](https://scienceparkstudygroup.github.io/rna-seq-lesson/setup.html).
> 2. Read the workshop [Code of Conduct](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html) to make sure this workshop stays welcoming for everybody.
> 3. Get comfortable: if you're not in a physical workshop, be set up with two screens if possible. You will be following along in RStudio on your own computer while also following this tutorial on your own.
> More instructions are available on the workshop website in the **Setup** section.
{: .prereq}

{% include links.md %}

## Citation
If you make use of this material in some way (teaching, vocational training, research), please cite us:
"Bliek Tijs, Frans van der Kloet and Marc Galland" (eds): "RNA-seq lesson."  Version 2020.04.
https://github.com/ScienceParkStudyGroup/rnaseq-lesson

## Credits
This lesson is heavily based on teaching materials from the [Harvard Chan Bioinformatics Core (HBC) in-depth NGS data analysis course](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/). Materials have been adapted and some exercises created to comply with the [Carpentries Foundation teaching requirements](https://carpentries.github.io/instructor-training/).


----



