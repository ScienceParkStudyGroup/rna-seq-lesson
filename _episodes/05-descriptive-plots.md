---
title: "Exploration of RNA-seq count results"
teaching: 60	
exercises: 30 
questions:
- "How are gene expression levels distributed within a RNA-seq experiment?"
- "Why do I need to scale/normalize read counts?"
- "How do I know that my RNA-seq experiment has worked according to my experimental design?"
- "How informative is PCA and sample clustering for RNA-seq quality checks?"
objectives:
- "Be able to calculate size factors and normalize counts using `DESeq2`."
- "Be able to create PCA and sample clustering plots to explore RNA-seq count results."
- "Be able to interpret a PCA plot and discuss its relationship with the experimental design."
keypoints:
- "Several biaises including sequencing depth can result in analysis artifacts and must be corrected trough scaling/normalisation."
- "RNA-seq results in a multivariate output that can be explored through data reduction methods (e.g. PCA)."
- "Sample clustering and PCA should indicate whether the observed experimental variability can be explained by the experimental design."
---

# Table of contents
<!-- MarkdownTOC autolink="True" levels="1,2" -->

- [1. Introduction](#1-introduction)
- [2. Normalization](#2-normalization)
	- [2.1 Common normalization methods](#21-common-normalization-methods)
	- [2.2 RPKM/FPKM \(not recommended for between sample comparisons\)](#22-rpkmfpkm-not-recommended-for-between-sample-comparisons)
	- [2.3 DESeq2-normalized counts: Median of ratios method](#23-deseq2-normalized-counts-median-of-ratios-method)
- [3. DESeq2 count normalization](#3-deseq2-count-normalization)
	- [3.1 Data import](#31-data-import)
	- [3.2 Match the experimental design and counts data](#32-match-the-experimental-design-and-counts-data)
	- [3.3 Create the DESeqDataSet object](#33-create-the-deseqdataset-object)
	- [3.4 Generate normalized counts](#34-generate-normalized-counts)
	- [3.5 Comparison of raw and scaled counts](#35-comparison-of-raw-and-scaled-counts)
- [4. Sample clustering](#4-sample-clustering)
	- [4.1 Distance calculation](#41-distance-calculation)
	- [4.2 Hierarchical clustering](#42-hierarchical-clustering)
	- [4.3 Heatmap of sample-to-sample correlations](#43-heatmap-of-sample-to-sample-correlations)
- [5. Principal Component Analysis](#5-principal-component-analysis)
	- [5.1 Introduction to PCA](#51-introduction-to-pca)
	- [5.2 The Iris data set](#52-the-iris-data-set)
	- [5.3 PCA applied to RNA-seq](#53-pca-applied-to-rna-seq)
- [6. Bonus: home-made DESeq normalization function](#6-bonus-home-made-deseq-normalization-function)
- [References](#references)
	- [Useful links](#useful-links)
	- [Photo credits](#photo-credits)

<!-- /MarkdownTOC -->


<img src="../img/05-workflow-overview.png" width="500px">


# 1. Introduction
The ultimate goal of most RNA-seq experiments is to accurately quantify the different transcripts present in a biological sample of interest. Yet, due to technical and biological causes, RNA-seq is prone to several biases that can affect sample / condition comparisons and therefore result interpretation.

# 2. Normalization

To be able to perform some quality check of our RNA-seq experiment and relate it to the experimental design, we will:
1. Normalize (scale) the gene counts obtained upon completion of the bioinformatic workflow. 
2. Perform some 
, which is necessary to make accurate comparisons of gene expression between samples.


The counts of mapped reads for each gene is proportional to the expression of RNA ("interesting") in addition to many other factors ("uninteresting"). Normalization is the process of scaling raw count values to account for the "uninteresting" factors. In this way the expression levels are more comparable between and/or within samples.

The main factors often considered during normalization are:

 - **Sequencing depth:** Accounting for sequencing depth is necessary for comparison of gene expression between samples. In the example below, each gene appears to have doubled in expression in *Sample A* relative to *Sample B*, however this is a consequence of *Sample A* having double the sequencing depth. 

    <img src="../img/normalization_methods_depth.png" width="400">

  >***NOTE:** In the figure above, each pink and green rectangle represents a read aligned to a gene. Reads connected by dashed lines connect a read spanning an intron.*

 - **Gene length:** Accounting for gene length is necessary for comparing expression between different genes within the same sample. In the example, *Gene X* and *Gene Y* have similar levels of expression, but the number of reads mapped to *Gene X* would be many more than the number mapped to *Gene Y* because *Gene X* is longer.

    <img src="../img/normalization_methods_length.png" width="200">

 - **RNA composition:** A few highly differentially expressed genes between samples, differences in the number of genes expressed between samples, or presence of contamination can skew some types of normalization methods. Accounting for RNA composition is recommended for accurate comparison of expression between samples, and is particularly important when performing differential expression analyses [[1](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)]. 

  In the example, imagine the sequencing depths are similar between Sample A and Sample B, and every gene except for gene DE presents similar expression level between samples. The counts in Sample B would be greatly skewed by the DE gene, which takes up most of the counts. Other genes for Sample B would therefore appear to be less expressed than those same genes in Sample A. 

<img src="../img/normalization_methods_composition.png" width="400">
    
***While normalization is essential for differential expression analyses, it is also necessary for exploratory data analysis, visualization of data, and whenever you are exploring or comparing counts between or within samples.***

## 2.1 Common normalization methods

Several common normalization methods exist to account for these differences:

| Normalization method | Description | Accounted factors | Recommendations for use |
| ---- | ---- | ---- | ---- |
| **CPM** (counts per million) | counts scaled by total number of reads | sequencing depth | gene count comparisons between replicates of the same samplegroup; **NOT for within sample comparisons or DE analysis**  |
| **TPM** (transcripts per kilobase million) | counts per length of transcript (kb) per million reads mapped | sequencing depth and gene length | gene count comparisons within a sample or between samples of the same sample group; **NOT for DE analysis** |
| **RPKM/FPKM** (reads/fragments per kilobase of exon per million reads/fragments mapped) | similar to TPM | sequencing depth and gene length | gene count comparisons between genes within a sample; **NOT for between sample comparisons or DE analysis** |
| DESeq2's **median of ratios** [[1](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)] | counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene | sequencing depth and RNA composition | gene count comparisons between samples and for **DE analysis**; **NOT for within sample comparisons** |
| EdgeR's **trimmed mean of M values (TMM)** [[2](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)] | uses a weighted trimmed mean of the log expression ratios between samples | sequencing depth, RNA composition, and gene length | gene count comparisons between and within samples and for **DE analysis** |

## 2.2 RPKM/FPKM (not recommended for between sample comparisons)

While TPM and RPKM/FPKM normalization methods both account for sequencing depth and gene length, RPKM/FPKM are not recommended. **The reason  is that the normalized count values output by the RPKM/FPKM method are not comparable between samples.** 

Using RPKM/FPKM normalization, the total number of RPKM/FPKM normalized counts for each sample will be different. Therefore, you cannot compare the normalized counts for each gene equally between samples. 

**RPKM-normalized counts table**

| gene   | sample_A | sample_B |
| ------ |:-----  : |:-------: |
| gene_1 | 5.5      | 5.5      |
| gene_2 | 73.4     | 21.8     |
| ...    | ...      | ...      |
|Total RPKM-normalized counts  | 1,000,000 | 1,500,000 |

For example, in the table above, Sample_A has a greater proportion of counts associated with gene_1  (5.5/1,000,000) than does sample_B (5.5/1,500,000) even though the RPKM count values are the same. Therefore, we cannot directly compare the counts for gene_1 (or any other gene) between sample_A and sample_B because the total number of normalized counts are different between samples. 

> *NOTE:* [This video by StatQuest](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) shows in more detail why TPM should be used in place of RPKM/FPKM if needing to normalize for sequencing depth and gene length.

## 2.3 DESeq2-normalized counts: Median of ratios method

Since tools for differential expression analysis are comparing the counts between sample groups for the same gene, gene length does not need to be accounted for by the tool. However, **sequencing depth** and **RNA composition** do need to be taken into account.

To normalize for sequencing depth and RNA composition, DESeq2 uses the median of ratios method. On the user-end there is only one step, but on the back-end there are multiple steps involved, as described below.

> **NOTE:**  The steps below describe in detail some of the steps performed by DESeq2 when you run a single function to get DE genes. Basically, for a typical RNA-seq analysis, **you would not run these steps individually** but rather make use of the `results()` function of `DESeq2`. 

**Step 1: creates a pseudo-reference sample (row-wise geometric mean)**

For each gene, a pseudo-reference sample is created that is equal to the geometric mean across all samples.

| gene | sample_A | sample_B | pseudo-reference sample  |
| ----- |:-----:|:-----:|:-----:|
| gene_1 | 1489 | 906 | sqrt(1489 * 906) = **1161.5** |
| gene_2 | 22 | 13 | sqrt(22 * 13) = **17.7** |
| ... | ... | ... | ... |

**Step 2: calculates ratio of each sample to the reference**

For every gene in a sample, the ratios (sample/ref) are calculated (as shown below). This is performed for each sample in the dataset. Since the majority of genes are not differentially expressed, the majority of genes in each sample should have similar ratios within the sample.

| gene | sample_A | sample_B | pseudo-reference sample  | ratio of sample_A/ref | ratio of sample_B/ref |
| ----- |:-----:|:-----:|:-----:| :-----: | :-----: |
| gene_1 | 1489 | 906 | 1161.5 | 1489/1161.5 = **1.28** | 906/1161.5 = **0.78** |
| gene_2 | 22 | 13 | 16.9 | 22/16.9 = **1.30** | 13/16.9 = **0.77** |
| gene_3 | 793 | 410 | 570.2 | 793/570.2 = **1.39** | 410/570.2 = **0.72**
| gene_4 | 76 | 42 | 56.5 | 76/56.5 = **1.35** | 42/56.5 = **0.74**
| gene_5 | 521 | 1196 | 883.7 | 521/883.7 = **0.590** | 1196/883.7 = **1.35** |
| ... | ... | ... | ... |

**Step 3: calculate the normalization factor for each sample (size factor)**

The median value (column-wise for the above table) of all ratios for a given sample is taken as the normalization factor (size factor) for that sample, as calculated below. Notice that the differentially expressed genes should not affect the median value:

`normalization_factor_sample_A <- median(c(1.28, 1.3, 1.39, 1.35, 0.59))`

`normalization_factor_sample_B <- median(c(0.78, 0.77, 0.72, 0.74, 1.35))`

The figure below illustrates the median value for the distribution of all gene ratios for a single sample (frequency is on the y-axis).

<img src="../img/deseq_median_of_ratios.png" width="400">

The median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample (large outlier genes will not represent the median ratio values). **This method is robust to imbalance in up-/down-regulation and large numbers of differentially expressed genes.**

> Usually these size factors are around 1, if you see large variations between samples it is important to take note since it might indicate the presence of extreme outliers.

**Step 4: calculate the normalized count values using the normalization factor**

This is performed by dividing each raw count value in a given sample by that sample's normalization factor to generate normalized count values. This is performed for all count values (every gene in every sample). For example, if the median ratio for sample_A was 1.3 and the median ratio for sample_B was 0.77, you could calculate normalized counts as follows:

sample_A median ratio = 1.3

sample_B median ratio = 0.77

**Raw Counts**

| gene | sample_A | sample_B |
| ----- |:-----:|:-----:|
| gene_1 | 1489 | 906 |
| gene_2 | 22 | 13 |
| ... | ... | ... |

**Normalized Counts**

| gene | sample_A | sample_B |
| ----- |:-----:|:-----:|
| gene_1 | 1489 / 1.3 = **1145.39** | 906 / 0.77 = **1176.62** |
| gene_2 | 22 / 1.3 = **16.92** | 13 / 0.77 = **16.88** |
| ... | ... | ... |

> ## Important note
> Please note that normalized count values are not integers anymore but rather decimal numbers. This is a good way to rapidly check whether you are dealing with normalised/scaled data and not raw gene counts. 
{: .callout}

***

# 3. DESeq2 count normalization 

## 3.1 Data import

~~~
# Data import 
counts <- read.delim("counts.txt", header = T, stringsAsFactors = F)
genes <- counts[,1]
counts <- counts[,-1]
row.names(counts) <- genes
xp_design <- read.delim("experimental_design_modified.txt", header = T, stringsAsFactors = F, colClasses = rep("character",4))

# change col names
colnames(xp_design) <- c("sample", "seed", "infected", "dpi")

# reorder counts columns according to the experimental design file
counts <- counts[, xp_design$sample]
~~~
{: .language-r}

This is how the first five rows/columns of the `counts` dataframe look like:
~~~
counts[1:5,1:5]
~~~
{: .language-r}

You can see that numbers in the matrix are not decimal numbers but rather integers. 
~~~
          ERR1406259 ERR1406271 ERR1406282 ERR1406294 ERR1406305
AT1G01010         59         79         61         71        123
AT1G01020        365        365        441        534        648
AT1G03987          8         14         14          4         20
AT1G01030        111        189        124        161        220
AT1G03993        131        155        183        236        250
~~~
{: .output}


Now that we know the theory of count normalization, we will normalize the counts using `DESeq2`. This requires a few steps:

1. Ensure the row names of the experimental design dataframe are present and in the same order as the column names of the counts dataframe.
2. Create a `DESeqDataSet` object.
3. Generate the normalized counts.

## 3.2 Match the experimental design and counts data

We should always make sure that we have sample names that match between the two files, and that the samples are in the right order. DESeq2 will output an error if this is not the case.

```r
## Check that sample names match in both files
all(colnames(counts) %in% xp_design$sample)
all(colnames(counts) == xp_design$sample)
```

If your data did not match, you could use the `match()` function to rearrange them to be matching.

## 3.3 Create the DESeqDataSet object

Bioconductor software packages often define and use a custom class within R for storing data (input data, intermediate data and also results). These custom data structures are similar to `lists` in that they can contain multiple different data types/structures within them. But, unlike lists they have pre-specified `data slots`, which hold specific types/classes of data. The data stored in these pre-specified slots can be accessed by using specific package-defined functions.

Let's start by creating the `DESeqDataSet` object and then we can talk a bit more about what is stored inside it. To create the object we will need the **count matrix** and the **experimental design** table as inputs. We will also need to specify a **design formula**. The design formula specifies the column(s) in the experimental design table and how they should be used in the analysis. For this tutorial, we only have one column of interest, that is `~ infected`. This column has two factor levels i.e. `mock` (not inoculated) and `Pseudomonas syringae_ DC3000` (inoculated). This column tells DESeq2 that for each gene we want to evaluate gene expression change with respect to these two different factor levels.


~~~
suppressPackageStartupMessages(library(DESeq2)) # to load DESeq2 and suppress the long startup message

# Creation of the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = xp_design, 
                              design = ~ 1) # the model does not matter here 
~~~
{: .language-r}

We now have a `DESeqDataSet` object that contains both count data and experimental metadata that is the relationship between samples and their combination of experimental factors. 

You can inspect this object by typing its name in your R console.
~~~
dds
~~~
{: .language-r}

~~~
class: DESeqDataSet 
dim: 33768 48 
metadata(1): version
assays(1): counts
rownames(33768): AT1G01010 AT1G01020 ... ATMG01400 ATMG01410
rowData names(0):
colnames(48): ERR1406259 ERR1406271 ... ERR1406285 ERR1406286
colData names(4): sample seed infected dpi
~~~
{: .output}

## 3.4 Generate normalized counts

The next step is to normalize the count data in order to be able to make fair gene comparisons between samples.


To perform the **median of ratios method** of normalization, DESeq2 has a single `estimateSizeFactors()` function that will generate size factors for us. We will use the function in the example below, but **in a typical RNA-seq analysis this step is automatically performed by the `DESeq()` function**, which we will see later. 

~~~
dds <- estimateSizeFactors(dds)
~~~
{: .language-r}

By assigning the results back to the `dds` object we are filling in the slots of the `DESeqDataSet` object with the appropriate information. We can take a look at the normalization factor applied to each sample using:

~~~
sizeFactors(dds)
~~~
{: .language-r}

We can plot these size factors to see how much they differ between samples. 
~~~
library(tidyverse)

# create a dplyr tibble
size_factors_df <- tibble(
  sample = names(sizeFactors(dds)), 
  correction_factor = sizeFactors(dds)
  )

# line plot to connect the different size factor values
p <- ggplot(size_factors_df, aes(x = sample, y = correction_factor, group = 1)) +
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(limits = c(0.5,2))

# to display the plot
p
~~~
{: .language-r}

This plot indicates that size factors are all between \~0.70 and \~1.8 so relatively close to each other. 

<img src="../img/05-size-factors.png" width="600px" alt="size factor plots" >

You can use DESeq-specific functions to access the different slots and retrieve information, if you wish. For example, suppose we wanted the original count matrix we would use `counts()`. For instance, to retrieve the normalized counts matrix from `dds`, we use the `counts()` function and add the argument `normalized=TRUE`.
~~~
# extract the normalised counts
counts_normalised = counts(dds, normalized = TRUE)
~~~
{: .language-r}

Let's take a peek at the first five rows and columns.
~~~~
counts_normalised[1:5,1:5]
~~~~
{: .language-r}


~~~
          ERR1406259 ERR1406271 ERR1406282 ERR1406294 ERR1406305
AT1G01010  72.152831   96.30085   66.47355  61.669649   88.51520
AT1G01020 446.369208  444.93431  480.57105 463.825247  466.32399
AT1G03987   9.783435   17.06597   15.25622   3.474346   14.39272
AT1G01030 135.745156  230.39064  135.12655 139.842443  158.31987
AT1G03993 160.203743  188.94471  199.42064 204.986439  179.90895
~~~
{: .output}
You now see that integers have become decimal numbers. All good!


> ## Note 
> `DESeq2` doesn't actually use normalized counts to compute differentially expressed genes. Rather, it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that peform differential expression analysis which use the negative binomial model.
{: .callout}

## 3.5 Comparison of raw and scaled counts

If you compare the raw and scaled gene count distributions, you should see that the medians (black bars) are more comparable after scaling. 

Here, only 6 samples out of the 48 samples are shown for clarity.    

<img src="../img/05-raw-scaled-counts.png" width="1000px" alt="Raw versus scaled counts">



# 4. Sample clustering
To assess if samples from the same condition are grouped together, we are going to perform a clustering analysis. It has two main steps:
- Calculate the distance between samples based on their normalized gene counts.
- Perform a hierarchical cluster analysis using `hclust`.
- Convert the result to a dendrogram and plot the resulting sample tree.  

## 4.1 Distance calculation
The Euclidean distance will be calculated between each sample based on the 
There are many ways to define a distance between two points.


Simple example first!
~~~
# sample 1 and 2
sample_a = c(1,2,3)
sample_b = c(1,2,3)

# Source: https://stackoverflow.com/questions/5559384/euclidean-distance-of-two-vectors
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
euc.dist(sample_a, sample_b)
~~~
{: .language-r}

The distance is 0 because the two samples have the same coordinates in their 3-dimensional space (x = 1, y = 1, z = 1).

~~~
sample_a = c(1,2,6)
sample_b = c(1,2,3)
euc.dist(sample_a, sample_b) # equal to 3

# if you do it manually
sqrt(sum(1 - 1)^2 + sum(2 - 2)^2 + sum(6 -3)^2)
~~~
{: .language-r}

This should give you the same result. 

Now, imagine you'd have to calculate the distance between each sample in a _N_ dimensional space (here N = 33,768 genes). Normally only Neo can see through the Matrix but maybe you do too...

<img src="../img/05-neo-matrix.jpg" width="600px">

Let's use R instead. 

We will work on the complete `counts_normalised` matrix as we only have to compute the distance between 48 samples. If you would do this for the \~33000 genes, it would take much much more time and you will have to stop your R session. 

~~~
# calculate the sample Euclidean distances 
distances_between_samples <- dist(
	t(counts_normalised),               # notice the t() for transpose. Otherwise you get the distance between genes
	method = "euclidean")  
as.matrix(distances_between_samples)[1:5,1:5]
~~~
{: .language-r}

~~~
           ERR1406259 ERR1406271 ERR1406282 ERR1406294 ERR1406305
ERR1406259        0.0   207590.4   167892.5   543575.1   823124.6
ERR1406271   207590.4        0.0   192718.5   417620.6   677362.6
ERR1406282   167892.5   192718.5        0.0   413844.6   692401.0
ERR1406294   543575.1   417620.6   413844.6        0.0   333604.2
ERR1406305   823124.6   677362.6   692401.0   333604.2        0.0
~~~
{: .output}

> ## Question
> How can you check that you have indeed calculated a distance between samples? 
> > ## Solution
> > You can notice that the distance between identical a sample and itself (e.g. ERR1406259 and ERR1406259) is equal to 0. 
> {: .solution}
{: .challenge}

## 4.2 Hierarchical clustering
This is the step where you will define your clusters. We will use the distance matrix computed before. 

There are different methods to perform clustering and they go way beyond this course but here is a note on the different methods ([see the original source](https://uc-r.github.io/hc_clustering#algorithms)).

> ## Methods
> **Maximum or complete linkage clustering:** It computes all pairwise dissimilarities between the elements in cluster 1 and the elements in cluster 2, and considers the largest value (i.e., maximum value) of these dissimilarities as the distance between the two clusters. It tends to produce more compact clusters.  
>
> **Minimum or single linkage clustering:** It computes all pairwise dissimilarities between the elements in cluster 1 and the elements in cluster 2, and considers the smallest of these dissimilarities as a linkage criterion. It tends to produce long, “loose” clusters.  
>
> **Mean or average linkage clustering:** It computes all pairwise dissimilarities between the elements in cluster 1 and the elements in cluster 2, and considers the average of these dissimilarities as the distance between the two clusters.  
>
> **Centroid linkage clustering:** It computes the dissimilarity between the centroid for cluster 1 (a mean vector of length p variables) and the centroid for cluster 2.  
>
> **Ward’s minimum variance method:** It minimizes the total within-cluster variance. At each step the pair of clusters with minimum between-cluster distance are merged.
{: .callout}

Here, we are going to use the Ward's clustering method. 
~~~
clustering_of_samples <- hclust(distances_between_samples, method = "ward.D2")
~~~
{: .language-r}

Finally, let's plot your results. 
~~~
plot(clustering_of_samples)
~~~
{: .language-r}

<img src="../img/05-dendro.png" width="600px">

## 4.3 Heatmap of sample-to-sample correlations

We can also build a more pretty figure where 
~~~
# to build a sample to 
xp_design_for_heatmap <- read.delim("experimental_design_modified.txt",
                                    header = TRUE,
                                    sep = "\t",
                                    colClasses = rep("factor",4))
                                    
    

row.names(xp_design_for_heatmap) <- xp_design_for_heatmap$sample
xp_design_for_heatmap$sample <- NULL

anno_info_colors = list(
  seed = c(MgCl2 = "#d8b365", 
           Methylobacterium_extorquens_PA1 = "#f5f5f5", 
           Sphingomonas_melonis_Fr1 = "#5ab4ac"),
  infected = c(mock = "gray", 
               Pseudomonas_syringae_DC3000 = "black"),
  dpi = c("2" = "dodgerblue", 
          "7" = "dodgerblue4")
)


correlation_between_samples <- cor(counts_normalised)

my_custom_breaks <- quantile(
  correlation_between_samples, # correlations
  seq(0,1,0.1))


# load library and plot heatmap
library(pheatmap)
library(RColorBrewer)

pheatmap(mat = correlation_between_samples, 
         color = rev(brewer.pal(n = length(my_custom_breaks) - 1,"RdYlBu")), 
         breaks = my_custom_breaks,
         annotation_row = xp_design_for_heatmap,
         annotation_col = xp_design_for_heatmap,
         annotation_colors = anno_info_colors
         )
~~~
{: .language-r}

<img src="../img/05-pheatmap.png" height="600px">

While the heatmap is full of colors and pleasant to display (arguable perhaps), it is not easy to distinguish between conditions. A Principal Component Analysis will help us to decide whether the experimental design is reflected in the count results.  

<br>
<br>

# 5. Principal Component Analysis

## 5.1 Introduction to PCA
In (bio)chemical analysis the data matrices can be very large. An infrared spectrum (800 wavelengths) for 50 samples for example would give a data matrix of size 40,000 (50x800) numbers.  A genomic data (e.g. 20,000 genes) for 100 patients would lead to a huge data matrix of (100x20,000) = 2,000,000 numbers. 

These matrices are so large that we need convenient ways to extract the important information from these large data matrices. 

Using principal component analysis (PCA) 

- The data is reduced to smaller matrices so they can more easily be examined, plotted and interpreted.
- The most important factors are extracted (principal components). These factors describe the multivariate (more than one variable) interactions between the measured variables.
- The samples can be classified to identify compound spectra, determine biomarkers etc.



To have an idea of how PCA works it should be noted that if there is a mutual relationship between two or more measurements (e.g. samples) they are correlated. These correlations can be strong (e.g. mass of object and weight on earth's surface) or weak (e.g. capabilities in sport and month of birth). In the example below there is a strong linear relationship between height and age for young children (*Moore, D.S. and McCabe G.P., Introduction to the Practice of Statistics (1989)*).

<img src="../img/age_height.png" alt="age_height" style="zoom:80%;" />

To explain PCA we use the example above and **project** the cases (observations) on a new coordinate system defined by principal component 1 (PC1) and principal component 2 (PC2) .

<img src="../img/age_height_to_scores.png"/>

In this new coordinate system, PC1 explains 99.77% of the total variation of the original data set while PC2 only explains 0.23%. Consequently, only variable (PC1) is sufficient to describe the whole data set which would support the conclusion that there is only 1 underlying factor, in this case age.

If the original data has more than two variables (e.g. n), which usually is the case, the projection would be in the n-dimensional space. Consequently more than two principal components can be calculated. By using an optimal projection, the principal components describe the maximum variance and are calculated in order of importance e.g.

![image-20200416141609987](../img/pc_exp_var_tbl.png)


## 5.2 The Iris data set

​The ability of PCA to capture as much variance as possible in the main principal components enables us to  to visualize (and summarize) the relation between objects (and variables) from a multi-dimensional space to a two dimensional plot.  



We can use the Fisher's famous Iris flower dataset from 1936 that describes three species of Iris (_Iris setosa_, _Iris virginica_ and _Iris versicolor_). The data set consists of 50 samples from each of these three species of Iris. Four features were measured from each sample: the length and the width of the sepals and petals, in centimeters ([source Wikipedia](https://en.wikipedia.org/wiki/Iris_flower_data_set)). See some pictures below.

| _Iris setosa_                                     | _Iris versicolor_                                  | _Iris virginica_                                        |
|---------------------------------------------------|----------------------------------------------------|---------------------------------------------------------|
| <img src="../img/Iris_setosa.jpg" width="200px"/> | <img src="../img/Iris_versicolor.jpg" width="200px"/> | <img src="../img/Iris_virginica.jpg" width="200px" />|


In the example below there is clear separation of the three types of [flowers](https://en.wikipedia.org/wiki/Iris_flower_data_set) which points to shared interactions between the different variables per group. 


To load this data in R  first the appropriate dataset have to be installed and loaded.
~~~
library(datasets)
data(iris)
head(iris)
~~~
{: .language-r}

~~~
Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
~~~
{: .output}

The summary of the iris data set display the content of the data. In this case the data consists of 150 objects with 5 variables. The last variable Species is a factor variable that specifies the class of the object.

~~~
 Sepal.Length    Sepal.Width     Petal.Length    Petal.Width          Species  
 Min.   :4.300   Min.   :2.000   Min.   :1.000   Min.   :0.100   setosa    :50  
 1st Qu.:5.100   1st Qu.:2.800   1st Qu.:1.600   1st Qu.:0.300   versicolor:50  
 Median :5.800   Median :3.000   Median :4.350   Median :1.300   virginica :50  
 Mean   :5.843   Mean   :3.057   Mean   :3.758   Mean   :1.199                  
 3rd Qu.:6.400   3rd Qu.:3.300   3rd Qu.:5.100   3rd Qu.:1.800                  
 Max.   :7.900   Max.   :4.400   Max.   :6.900   Max.   :2.500 
~~~
{: .output}


For convenience we use a very rudimentary (own) implementation implementation of PCA. Copy-paste this code into R to load this function into your environment and use it later on.

~~~
# define a custom R function called "mypca()""
mypca <- function(x, center = TRUE, scale = TRUE){
  
  # This checks that the samples/individuals are in rows
  # If samples in columns, then the matrix is transposed
  if(nrow(x)>ncol(x))
  {
    x = t(x)
  }
  
  # remove constant variables
  constant_val = apply(x,2,'sd')
  x_reduced = x[,constant_val>0]
  
  # perform SVD
  SVD <- svd(scale(x_reduced,center = center, scale = scale))
  
  # create scores data frame
  scores <- as.data.frame(SVD$u %*% diag(SVD$d))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC", c(1:dim(scores)[2]))
  
  # create loadings data frams
  loadings <- data.frame(SVD$v)
  colnames(loadings) <- paste0("PC", c(1:dim(loadings)[2]))
  rownames(loadings) <- colnames(x_reduced)
  
  # create data frame for explained variances
  explained_var <- as.data.frame(round((SVD$d^2) / sum(SVD$d^2)*100, digits = 1))
  rownames(explained_var) <- paste0("PC", c(1:dim(loadings)[2]))
  colnames(explained_var) <- "exp_var"
  
  # return result
  return (list("scores" = scores, "loadings" = loadings, "explained_var" = explained_var))
}
~~~
{: .language-r}

Now we have everything in our R environment into place, we can actually perform the PCA analysis and create the plots.  
Since the four first principal components catch most if not all 

~~~
# perform the PCA analysis on only the first 4 variables (skip the Species variable)
pca <- mypca(iris[,1:4], scale = FALSE)
scores = as.data.frame(pca$scores[,1:2])
scores['Species'] = iris$Species
explained_var = pca$explained_var$exp_var
~~~
{: .language-r}

In order to have an idea of how effective the 'compression' or variable reduction of the PCA algorithm was on our data set, we make a so-called 'scree' plot in which the explained variance is expressed as a function of the number of principal components.

~~~
# add a convenient column number for the bar plot to display
dfev <- data.frame(PC = c(1,2,3,4), exp_var  = pca$explained_var)

# make the plot
scree_plot <- ggplot(dfev, aes(x = PC, y = exp_var)) +
       ylab('explained variance (%)') + 
       ggtitle('explained variance per component') + 
       geom_bar(stat = "identity")

# display it
scree_plot
~~~
{: .language-r}

<img src="../img/pca_iris_exp_var.png" width="600px" alt="pca_iris_exp_var PCA">

It is clear that the PCA really reduces our data to almost 1 variable (component) which explains over 90% of all the variance of the data. The second component only explains 5.3%. The third and fourth even less.

The whole idea behind the analysis is to visualize the high-dimensional data (e.g. a score plot) in a smaller dimensional space (e.g. 2D space). In R we can do that with the following lines of code

~~~
# plot the scores of the first 2 components
p <- ggplot(scores) + 
       geom_point(aes(x = PC1, y = PC2, shape = Species, col = Species)) + 
       xlab(paste0('PC1(',explained_var[1],'%)')) + 
       ylab(paste0('PC2(',explained_var[2],'%)')) + 
       ggtitle('PCA score plot')
p
~~~
{: .language-r}

<img src="../img/pca_iris_new.png" width="600px" alt="pca_iris_new">

From the score plot it is clear that the Setosa flowers are clearly different from the Versicolor/Virginica flowers. Versicolor and Virginica cannot be separated on PC1 and/or PC2. Looking at the PC1 vs PC3 however, the two groups can be separated better. It is very important to understand that even if a principal component explains a low amount of variance it still can contain interesting (biological) information. 

<img src="../img/pca_iris_1_3.png" width="600px" alt="pca_iris_1_3">

The scores are indicative of how the objects in the data set score in the new component space, correspondingly the loadings indicate how the variables score in the component space. The score plots above for example show a separation on PC1 between the 3 groups. If we would like to know which variables are important for this separation we can try to interpret our data.



~~~
library(reshape2) # to access the melt() function

# reformat the loading data
loadings <- melt(pca$loadings)
# rename the columns
colnames(loadings)<-c("Component","Value")
# add the 'original' variable names
loadings['Variable']=as.factor(rep(colnames(iris)[-5],4))
# plot the loading values per components
p <- ggplot(loadings,aes(x=Variable,y=Value)) +  geom_bar(stat='identity') + facet_wrap(~Component)
p
~~~
{: .language-r}


<img src="../img/pca_iris_loadings.png" width="600px" alt="Iris loadings">


From the loading plot for PC1 it is clear that `Petal.Length` is the most important factor while `Sepal.Width` is the one which is least important. Because PC1 explains 92.5 % of the total variance we can conclude that `Petal.Length` is the most important factor that separates the three groups of flowers.   
There are many more things to learn on PCA (e.g. scaling, filtering) but that is out of the scope of these exercises. 


## 5.3 PCA applied to RNA-seq

In the context of an RNA-seq experiment, it can be used to visualize the differences (distances) between samples and how it relates to the experimental design.

Let's plot the samples along the two first components.

This can be done with the `plotPCA()` function of the `DESeq2` package. First, we need to stabilise the variance across genes with different means using a variance stabilising transformation or `vst()` . If interested, you can check the [corresponding detailed section in the DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization).

First, let's load the library that we need. 
~~~
library(DESeq2, quietly = TRUE)
library(ggplot2)
~~~
{: .language-r}

Let's make a first version of a PCA plot using the `plotPCA()` function of the `DESeq2` library. We are going to plot our experimental design on several PCA plots. 

### 5.3.1 Infected versus mock
Let's first see how our _P. syringae_ infection condition is reflected at the PCA level.

~~~
# variance-stabilizing transformation
vst_dds <- vst(dds)

# PCA using the plotPCA function
pcaData <- plotPCA(vst_dds, intgroup = c("infected"), returnData = TRUE) 

percentVar <- round(100 * attr(pcaData, "percentVar"))

# customised PCA plot
ggplot(pcaData, aes(PC1, PC2, color = infected)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")
~~~
{: .language-r}

<img src="../img/05-pca-1.png" width="800px" alt="PCA with infected conditions overlayed">

This PCA plot clearly shows that the infected condition seems to be related to the first principal component (59% of the variance).

But we can also overlay different experimental conditions on this PCA plot. A PCA analysis let the data "self-organise".
Have a look at this new PCA plot with the seed conditions superimposed.  

~~~
# PCA using the plotPCA function
pcaData <- plotPCA(vst_dds, intgroup = c("seed"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# customised PCA plot
ggplot(pcaData, aes(PC1, PC2, color = seed)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")
~~~
{: .language-r}

<img src="../img/05-pca-2.png" width="800px" alt="PCA with seed inoculation conditions overlayed" >


You can see that the sample positions have not changed on the plot. But the colors have since we wanted to color samples by their seed inoculation level ("MgCl2", "Fr1" or "PA1").

We are now creating a final and complete version of the PCA with all factors included. It is quite a complex plot but you should find previous patterns back.
~~~
# PCA using the plotPCA function
pcaData <- plotPCA(vst_dds, intgroup = c("seed", "infected", "dpi"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

# customised PCA plot
ggplot(pcaData, aes(PC1, PC2, color = seed, shape = infected, size = dpi)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")
~~~
{: .language-r}

<img src="../img/05-pca-3.png" width="800px" alt="complete PCA" >

# 6. Bonus: home-made DESeq normalization function

Please consult the step-by-step R code to normalize the DESeq2 way [here](../median_of_ratios_manual_normalization/index.html).

# References

## Useful links
1. [Gabriel Martos cluster analysis](https://rpubs.com/gabrielmartos/ClusterAnalysis)
2. [Bradley Boehmke](https://uc-r.github.io/hc_clustering)
3. Love et al. (2014) [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. _Genome Biology_, 15:550](https://link.springer.com/article/10.1186/s13059-014-0550-8).
4. [Statquest](https://www.youtube.com/watch?v=UFB993xufUU)
5. [Harvard Bioinformatic Core Training program](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

## Photo credits
<a style="background-color:black;color:white;text-decoration:none;padding:4px 6px;font-family:-apple-system, BlinkMacSystemFont, &quot;San Francisco&quot;, &quot;Helvetica Neue&quot;, Helvetica, Ubuntu, Roboto, Noto, &quot;Segoe UI&quot;, Arial, sans-serif;font-size:12px;font-weight:bold;line-height:1.2;display:inline-block;border-radius:3px" href="https://unsplash.com/@markusspiske?utm_medium=referral&amp;utm_campaign=photographer-credit&amp;utm_content=creditBadge" target="_blank" rel="noopener noreferrer" title="Download free do whatever you want high-resolution photos from Markus Spiske"><span style="display:inline-block;padding:2px 3px"><svg xmlns="http://www.w3.org/2000/svg" style="height:12px;width:auto;position:relative;vertical-align:middle;top:-2px;fill:white" viewBox="0 0 32 32"><title>unsplash-logo</title><path d="M10 9V0h12v9H10zm12 5h10v18H0V14h10v9h12v-9z"></path></svg></span><span style="display:inline-block;padding:2px 3px">Markus Spiske</span></a>


