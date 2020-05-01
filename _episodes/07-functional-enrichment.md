---
title: "Going beyond a list of genes"
teaching: 30
exercises: 60 
questions:
- "What are factor levels and why is it important for different expression analysis?"
objectives:
- "Be able to perform a gene set enrichment analysis (GSEA) using a list of differentially expressed genes."
keypoints:
- "ORA can...."

---


# Table of Contents
<!-- MarkdownTOC autolink="True" levels="1,2,3" -->

- [1. Introduction](#1-introduction)
- [2. Annotating your DE genes](#2-annotating-your-de-genes)
  - [2.1 Load the table of differential genes](#21-load-the-table-of-differential-genes)
  - [2.2 Ensembl databases](#22-ensembl-databases)
  - [2.3 Querying Ensembl databases using biomartr](#23-querying-ensembl-databases-using-biomartr)
- [3. Gene set enrichment analysis and overrepresentation analysis](#3-gene-set-enrichment-analysis-and-overrepresentation-analysis)
  - [3.1 Over Representation Analysis \(ORA\)](#31-over-representation-analysis-ora)
  - [3.1 Using R](#31-using-r)
  - [3.3 Metascape](#33-metascape)
  - [3.2 AgriGO](#32-agrigo)
- [4. Data integration with metabolic pathways](#4-data-integration-with-metabolic-pathways)
  - [Using MapMan](#using-mapman)
  - [Using iPath](#using-ipath)
- [5. Looking for regulatory elements](#5-looking-for-regulatory-elements)
  - [5.1 Extracting the coordinqtes of genes](#51-extracting-the-coordinqtes-of-genes)
  - [5.2 Adding or substracting X nts](#52-adding-or-substracting-x-nts)
- [6. Other sources of information](#6-other-sources-of-information)
  - [6.1 ThaleMiner](#61-thaleminer)
  - [6.2 Expression atlas](#62-expression-atlas)
  - [6.3 BAR](#63-bar)
  - [6.4 CoExprViz](#64-coexprviz)
- [7. Useful links](#7-useful-links)
- [References](#references)

<!-- /MarkdownTOC -->



# 1. Introduction
You've finally managed to extract a list of differentially expressed genes from your comparison. Great job! But...now what? :question: :confused:

Why did you do the experiment in the first place? Probably because you had an hypothesis or you were looking to open new leads. 

In this tutorial, we are looking for Arabidopsis leaf genes that are induced or repressed upon inoculation by _Pseudomonas syringae_ DC3000 after 7 days.
One important goal is to gain a higher view and not only deal with individual genes but understand which pathways are involved in the response. 

Once we obtain a list of genes, we have multiple analysis to perform to go beyond a simple list of genes:
- a Gene Set Enrichment Analysis (GSEA)

- Co-expression analysis: can you identify genes working together upon your experimental treatment?
- Tissue-specificity 

# 2. Annotating your DE genes

We are going to use two fantastic resources: the superb [Ensembl](https://www.ensembl.org) databases and the marvelous [biomartr package](https://docs.ropensci.org/biomartr).
Together, they will automate a lot of tedious and tiring steps when you want to retrieve gene annotations, sequences, etc.


## 2.1 Load the table of differential genes

~~~
diff_genes <- read_delim(file = "03.RNA-seq/differential_genes.tsv", delim = "\t")
~~~
{: .language-r}

All what we know about the differential genes are their locus identifier. Not much....
We are missing functional information which we will add. 

## 2.2 Ensembl databases 
Gene information can be automatically queried directly from R to access the [Ensembl](https://www.ensembl.org) databases. 
Ensembl gathers a tremendous amount of genomic information which can be accessed through a web browser or programmatically. 

> The Ensembl (https://www.ensembl.org) is a system for generating and distributing genome annotation such as genes, variation, regulation and comparative genomics across the vertebrate subphylum and key model organisms. The Ensembl annotation pipeline is capable of integrating experimental and reference data from multiple providers into a single integrated resource. Here, we present 94 newly annotated and re-annotated genomes, bringing the total number of genomes offered by Ensembl to 227. 

<img src="../img/07-ensembl.png" alt="Portal of Ensembl for Arabidopsis" width="1000px">

## 2.3 Querying Ensembl databases using biomartr

What purpose serves `biomartr`? From the documentation:
> The first step, however, of any genome based study is to retrieve genomes and their annotation from databases. To automate the retrieval process of this information on a meta-genomic scale, the biomartr package provides interface functions for genomic sequence retrieval and functional annotation retrieval. The major aim of biomartr is to facilitate computational reproducibility and large-scale handling of genomic data for (meta-)genomic analyses. In addition, biomartr aims to address the genome version crisis. With biomartr users can now control and be informed about the genome versions they retrieve automatically. Many large scale genomics studies lack this information and thus, reproducibility and data interpretation become nearly impossible when documentation of genome version information gets neglected.

What is available for _Arabidopsis thaliana_ in Ensembl?
~~~
biomartr::organismBM(organism = "Arabidopsis thaliana")
~~~
{: .language-r}

~~~
organism_name description                                                  mart         dataset       version
  <chr>         <chr>                                                        <chr>        <chr>         <chr>  
1 athaliana     Arabidopsis thaliana genes (TAIR10)                          plants_mart  athaliana_eg… TAIR10 
2 athaliana     Arabidopsis thaliana Short Variants (SNPs and indels exclud… plants_vari… athaliana_eg… TAIR10 
~~~
{: .output}

This indicates that we can get a dataset called `athaliana_eg_gene` of the genome annotation version `TAIR10` from the `plant_mart` mart. 

Let's see how many different information fields we can retrieve from the `arabidopsis_eg_gene` dataset.
~~~
arabido_attributes = 
  biomartr::organismAttributes("Arabidopsis thaliana") %>% 
  filter(dataset == "athaliana_eg_gene")
arabido_attributes
~~~
{: .language-r}


~~~
# A tibble: 2,574 x 4
   name                  description              dataset           mart       
   <chr>                 <chr>                    <chr>             <chr>      
 1 ensembl_gene_id       Gene stable ID           athaliana_eg_gene plants_mart
 2 ensembl_transcript_id Transcript stable ID     athaliana_eg_gene plants_mart
 3 ensembl_peptide_id    Protein stable ID        athaliana_eg_gene plants_mart
 4 ensembl_exon_id       Exon stable ID           athaliana_eg_gene plants_mart
 5 description           Gene description         athaliana_eg_gene plants_mart
 6 chromosome_name       Chromosome/scaffold name athaliana_eg_gene plants_mart
 7 start_position        Gene start (bp)          athaliana_eg_gene plants_mart
 8 end_position          Gene end (bp)            athaliana_eg_gene plants_mart
 9 strand                Strand                   athaliana_eg_gene plants_mart
10 band                  Karyotype band           athaliana_eg_gene plants_mart
# … with 2,564 more rows
~~~
{: .output}

:scream: There is quite some information in there! We should be able to get what we want!

~~~
attributes_to_retrieve = c("tair_symbol", "uniprotswissprot", "kegg_enzyme")

result_BM <- biomartr::biomart( genes      = diff_genes$genes,                  # genes were retrieved using biomartr::getGenome()
                                mart       = "plants_mart",                     # marts were selected with biomartr::getMarts()
                                dataset    = "athaliana_eg_gene",               # datasets were selected with biomartr::getDatasets()
                                attributes = attributes_to_retrieve,            # attributes were selected with biomartr::getAttributes()
                                filters = 	"ensembl_gene_id" )# query key
head(result_BM)  
~~~
{: .language-r}

# 3. Gene set enrichment analysis and overrepresentation analysis

## 3.1 Over Representation Analysis (ORA)

Over Representation Analysis is searching for biological functions or pathways that are enriched in a list obtained through experimental studies compared to the complete list of functions/pathways.  


$$p = 1 - {\sum_{i=0}^{k-1} {M \choose i}{N - M \choose n - i} \over {N \choose n}}$$  

In this formula: 
- **N** is the total number of genes in the background distribution.
- **M** is the number of genes within that distribution that are annotated (either directly or indirectly) to the gene set of interest.
- **n** is the size of the list of genes of interest (the size of your "drawing").
- **k** and k is the number of genes within that list which are annotated to the gene set. The background distribution by default is all the genes that have annotation. P-values should be adjusted for multiple comparison.

Do you remember your math classes from high school? Now's the time to get them to work again!

Binomial coefficient is defined as $${n \choose k}$$ and is equal to $$n! \over {k! (n-k)!}$$



## 3.1 Using R
Cluster profiler: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html




## 3.3 Metascape

## 3.2 AgriGO
Agrigo


# 4. Data integration with metabolic pathways

ClusterProfiler: https://yulab-smu.github.io/clusterProfiler-book/chapter6.html

## Using MapMan 
prep of the data table
visualisation on pathways

From [Schwacke et al., 2019](https://doi.org/10.1016/j.molp.2019.01.003):
> The MapMan framework was developed specifically for plants with the design goal to facilitate the visualization of omicsdata on plant pathways (Thimm et al., 2004). 
> MapMan uses a simple hierarchical tree structure of terms referred to as ‘‘bins,’’ which describe biological contexts/concepts.   
> Major biological processes (e.g., photosynthesis) are encompassed in top-level bins, and each child bin represents a more narrowly focused subprocess or component 
> within the context of the parent bin. Assignment of proteins to the lowest-level (i.e., leaf) bins was preferred in order to make the annotation as precise as possible, although assignment to abstract higher-level bins was supported.   
> Proteins were mostly assigned to a single bin, but for some proteins with functions in diverse biological processes it wasnecessary to correspondingly assign to multiple bins.

## Using iPath
KeggKOALA 

# 5. Looking for regulatory elements

## 5.1 Extracting the coordinqtes of genes

## 5.2 Adding or substracting X nts
For instance, 5000 nts 
If gene is on DNA strand + then substract 5000 nts  
If gene is on DNA strand - then add 5000 nts

Promoter retrieval using GenomicRanges
MEME for motif...

# 6. Other sources of information

## 6.1 ThaleMiner
https://bar.utoronto.ca/thalemine/begin.do

## 6.2 Expression atlas
https://www.ebi.ac.uk/gxa/home

## 6.3 BAR
http://www.bar.utoronto.ca/

## 6.4 CoExprViz
http://bioinformatics.psb.ugent.be/webtools/coexpr/

# 7. Useful links
- [BiomartR](https://docs.ropensci.org/biomartr/)
- [Arabidopsis.org (TAIR) list of data mining tools](https://www.arabidopsis.org/portals/expression/microarray/microarrayExpressionV2.jsp)
- [ResearchGate related question](https://www.researchgate.net/post/How_can_I_analyze_a_set_of_DEGs_differentially_expressed_genes_to_obtain_information_from_them)	


# References
* [The Cluster Profiler companion boo, a great place to start](https://yulab-smu.github.io/clusterProfiler-book/chapter2.html)
* Zhou et al. (2019). Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nat Commun 10, 1523 (2019). [link](https://doi.org/10.1038/s41467-019-09234-6)
* Yates et al. (2020) Ensembl 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D682–D688, [Link](https://doi.org/10.1093/nar/gkz966)
* Tian et al. (2017) agriGO v2.0: a GO analysis toolkit for the agricultural community. _Nucleic Acids Research_, Volume 45, Issue W1, Pages W122–W129.[Link](https://doi.org/10.1093/nar/gkx382) 
* MapMan: [MapMan4: A Refined Protein Classification and Annotation Framework Applicable to Multi-Omics Data Analysis. Schwacke et al. _Molecular Plant_, 12(6):879-892](https://doi.org/10.1016/j.molp.2019.01.003)
* Drost et al. (2017) Biomartr: genomic data retrieval with R. _Bioinformatics_ 33(8): 1216-1217. [doi:10.1093/bioinformatics/btw821](https://academic.oup.com/bioinformatics/article/33/8/1216/2931816.




