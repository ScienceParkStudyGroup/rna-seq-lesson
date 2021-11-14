---
title: "09 Transcriptomic and metabolomic data integration"
teaching: 45
exercises: 15 
questions:
- "How can I map differential genes to metabolic pathways?"
- "How do I retrieve KEGG identifiers given a list of gene identifiers?"
objectives:
- "Retrieve common metabolic identifiers such as KEGG KO using the Ensembl API queried using the `biomartr` package. "
- "Be able to use MapMan to generate pathway visualisations."
keypoints:
- "Transcriptomic data integration with metabolic pathways require to map gene identifiers to pathway  the correspondence "

---

# Table of Contents
<!-- MarkdownTOC autolink="True" levels="1,2" -->

- [1. Introduction](#1-introduction)
  - [1.2 Setup](#12-setup)
- [2. iPath webtool](#2-ipath-webtool)
  - [2.1 Tool presentation](#21-tool-presentation)
  - [2.2 Retrieving Uniprot identifiers](#22-retrieving-uniprot-identifiers)
  - [5.2 MapMan](#52-mapman)
- [6. Other data mining tools](#6-other-data-mining-tools)
  - [6.1 ThaleMiner](#61-thaleminer)
  - [6.2 Expression atlas](#62-expression-atlas)
  - [6.3 BAR](#63-bar)
  - [6.4 CoExprViz](#64-coexprviz)
- [7. Troubleshooting](#7-troubleshooting)
- [8. Going further](#8-going-further)
  - [8.1 Useful links](#81-useful-links)
  - [8.2. References](#82-references)

<!-- /MarkdownTOC -->

<img src="../img/07-workflow-overview.png" width="500px" alt="workflow overview">

# 1. Introduction



## 1.2 Setup

~~~

~~~
{: .language-r}


# 2. iPath webtool 

## 2.1 Tool presentation
So far, we have only been mostly looking at our transcriptomic results alone. Yet, KEGG has started to give us some insight into metabolism. Let's dig further with two tools: **iPath** and  **MapMan**.  

From the [Interactive Pathways Explorer v3 home page](https://pathways.embl.de/):
> Interactive Pathways Explorer (iPath) is a web-based tool for the visualization, analysis and customization of various pathway maps.  
iPath provides extensive map customization and data mapping capablities. Colors, width and opacity of any map element can be changed using various types of data (for example KEGG KOs, COGs or EC numbers). 

We will use this webtool to map our genes on metabolic pathways to provide some contextual information. 

<img src="../img/07-ipath-1.png" alt="home page of iPath 3 web interface" height="400px">

## 2.2 Retrieving Uniprot identifiers 

From our list of differential genes
We are going to generate a list of Uniprot identifiers since these are usable to visualise with iPath3 from our list of differential genes.
~~~
diff_arabidopsis_genes_annotated %>% 
  filter(uniprotswissprot != "") %>%                                       # to remove genes with no matching Uniprot entries
  unique() %>% 
  mutate(id_for_ipath = paste("UNIPROT",uniprotswissprot,sep = ":")) %>%   # to create an ID that iPath can use
  dplyr::select(id_for_ipath) %>%                                          # we keep only the relevant ID for further copy-pasting 
  write.table(., 
    file = "diff_genes_swissprot.tsv", 
    row.names = FALSE, 
    quote = FALSE)
~~~
{: .language-r}

If you open the `diff_genes_swissprot.tsv` file in a text editor, this is what you should see. 
~~~
id_for_ipath
UNIPROT:Q9MAN1
UNIPROT:Q5XEZ0
UNIPROT:O24457
UNIPROT:Q9MAM1
UNIPROT:Q9LNJ9
UNIPROT:Q06402
... more lines ....
~~~
{: .output}

We can copy-paste this list in the box to the right. 

<img src="../img/07-ipath-3.png" height="400px">

If you then click on the "Submit data" button, you should see the following reactions highlighted.  

<img src="../img/07-ipath-2.png" alt="first metabolic pathway map" height="400px">

Ok, there is a lot of information there and possibly too much. Let's filter out some genes based on their $$\log_{2}$$ fold change. 

We can calculate the median, 75th percentile and the 90th percentile of our fold changes. 
~~~
diff_genes %>% 
  filter(log2FoldChange > 0) %>% 
  with(.,quantile(log2FoldChange, c(0.5,0.75,0.9)))
~~~
{: .language-r}

~~~
50%       75%       90% 
0.8695469 1.6819535 3.2561618 
~~~
{: .output}
Since 75% of the genes with a positive $$\log_{2}$$ fold change have a value lower than 1.68, we are going to keep the genes with a higher fold change. 

~~~
diff_genes_filtered = 
  diff_genes %>% 
  filter(log2FoldChange > quantile(log2FoldChange, 0.75)) 

# we query Ensembl again to retrieve the attributes
# attributes_to_retrieve = c("tair_symbol", "uniprotswissprot","entrezgene_id")
diff_arabidopsis_genes_annotated_2 <- biomartr::biomart(genes = diff_genes_filtered$genes,
                                                      mart       = "plants_mart",                 
                                                      dataset    = "athaliana_eg_gene",           
                                                      attributes = attributes_to_retrieve,        
                                                      filters =     "ensembl_gene_id" )  

diff_arabidopsis_genes_annotated_2 %>% 
  filter(uniprotswissprot != "") %>% 
  unique() %>% 
  mutate(id_for_ipath = paste("UNIPROT",uniprotswissprot,sep = ":")) %>% 
  dplyr::select(id_for_ipath) %>% 
  write.table(., file = "diff_genes_swissprot_2.tsv", row.names = FALSE, quote = FALSE)
~~~
{: .language-r}

In your `diff_genes_swissprot_2.tsv` file, you should have less entries now which will make the map slightly clearer. 

<img src="../img/07-ipath-4.png" alt="simplified ipath metabolic map" height="400px">


## 5.2 MapMan
 
MapMan is a standalone software that you can install on your laptop to create insightful metabolic visualisation of transcriptomic, proteomic and/or metabolic data. 

From [Schwacke et al., 2019](https://doi.org/10.1016/j.molp.2019.01.003):
> The MapMan framework was developed specifically for plants with the design goal to facilitate the visualization of omicsdata on plant pathways (Thimm et al., 2004). 
> MapMan uses a simple hierarchical tree structure of terms referred to as ‘‘bins,’’ which describe biological contexts/concepts.   
> Major biological processes (e.g., photosynthesis) are encompassed in top-level bins, and each child bin represents a more narrowly focused subprocess or component 
> within the context of the parent bin. Assignment of proteins to the lowest-level (i.e., leaf) bins was preferred in order to make the annotation as precise as possible, although assignment to abstract higher-level bins was supported.   
> Proteins were mostly assigned to a single bin, but for some proteins with functions in diverse biological processes it wasnecessary to correspondingly assign to multiple bins.

### 5.2.1 Requirements
You will need to create an account first and download the proper MapMan installation for your OS [here](https://mapman.gabipd.org/mapman-download). I'd recommend to download the latest stable version (not the "release candidate").

You will then have to download the so-called "mapping file" which is an ontology relating genes to functional bins. You can also generate your own classification for your species of interest with the [Mercator webtool](https://www.plabipd.de/portal/mercator4). You'll need to supply a fasta file containing your transcripts or your protein of interest. 

The latest list of available mappings is [here](https://mapman.gabipd.org/mapmanstore?p_p_id=MapManDataDownload_WAR_MapManDataDownloadportlet_INSTANCE_4Yx5&p_p_lifecycle=0&p_p_state=normal&p_p_mode=view&p_p_col_id=column-1&p_p_col_pos=1&p_p_col_count=2). 


### 5.2.2 Pathway visualisation
One of the most used visualisation is the global overview of the metabolism. Up-regulated genes are shown in red and genes down-regulated in blue. 

<img src="../img/07-mapman-1.jpeg" alt="mapman visualisation" height="400px">

# 6. Other data mining tools

## 6.1 ThaleMiner
[https://bar.utoronto.ca/thalemine/begin.do](https://bar.utoronto.ca/thalemine/begin.do)

## 6.2 Expression atlas
[https://www.ebi.ac.uk/gxa/home](https://www.ebi.ac.uk/gxa/home)

## 6.3 BAR
[http://www.bar.utoronto.ca/](http://www.bar.utoronto.ca/)

## 6.4 CoExprViz
[http://bioinformatics.psb.ugent.be/webtools/coexpr/](http://bioinformatics.psb.ugent.be/webtools/coexpr/)

# 7. Troubleshooting
If biomart refuses to query Ensembl again, run this command:
~~~
biomaRt::biomartCacheClear() # to solve a known bug https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/335
~~~
{: .language-r}

This will clean the cache memory and allow to perform the Ensembl query again.  
# 8. Going further 

## 8.1 Useful links
- [BiomartR](https://docs.ropensci.org/biomartr/)
- [Arabidopsis.org (TAIR) list of data mining tools](https://www.arabidopsis.org/portals/expression/microarray/microarrayExpressionV2.jsp)
- [ResearchGate related question](https://www.researchgate.net/post/How_can_I_analyze_a_set_of_DEGs_differentially_expressed_genes_to_obtain_information_from_them)	

## 8.2. References
* [The Cluster Profiler companion book, a great place to start](https://yulab-smu.github.io/clusterProfiler-book/chapter2.html)
* Zhou et al. (2019). Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nat Commun 10, 1523 (2019). [link](https://doi.org/10.1038/s41467-019-09234-6)
* Yates et al. (2020) Ensembl 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D682–D688, [Link](https://doi.org/10.1093/nar/gkz966)
* Tian et al. (2017) agriGO v2.0: a GO analysis toolkit for the agricultural community. _Nucleic Acids Research_, Volume 45, Issue W1, Pages W122–W129.[Link](https://doi.org/10.1093/nar/gkx382) 
* MapMan: [MapMan4: A Refined Protein Classification and Annotation Framework Applicable to Multi-Omics Data Analysis. Schwacke et al. _Molecular Plant_, 12(6):879-892](https://doi.org/10.1016/j.molp.2019.01.003)
* Drost et al. (2017) Biomartr: genomic data retrieval with R. _Bioinformatics_ 33(8): 1216-1217. [doi:10.1093/bioinformatics/btw821](https://academic.oup.com/bioinformatics/article/33/8/1216/2931816).
* Darzi et al. (2018) iPath3.0: interactive pathways explorer v3. _Nucleic Acids Research_, Volume 46, Issue W1, 2 July 2018, Pages W510–W513, [link](https://doi.org/10.1093/nar/gky299)


