---
title: "Differential expression analysis"
teaching: 30
exercises: 60 
questions:
- "Q1"
- "Q2"
- "Q3"
- "Q4"
objectives:
- "O1"
- "O2"
- "O3"
keypoints:
- "KP1"
- "KP2"
- "KP3"

---



## Table of Contents
1. [Introduction](#introduction)
2. [Creating the DESeqDataSet object](#creating-the-DESeqDataSet-object)
3. [Different expression analysis](#creating-the-DESeqDataSet-object)



## Introduction

In this episode, we will see how to perform a simple one-condition experimental comparison. We will compare the transcriptome of Arabidopsis in response to infection by the leaf pathogenic bacteria _Pseudomonas syringae DC3000_ after 7 days (7 dpi). 


## Creating the DESeqDataSet object
We are first going to import the count file and the experimental design file required to build
a DESeqDataSet object from input tables. 

~~~
# import counts and experimental design 
counts <- read.delim("counts.txt", header = T, stringsAsFactors = F)
genes <- counts[,1]
counts <- counts[,-1] # to eliminate the gene columns that won't match sample names in the design file.

xp_design <- read.delim("experimental_design_modified.txt", header = T, stringsAsFactors = F, colClasses = rep("character",4))

# change col names
colnames(xp_design) <- c("sample", "growth", "infected", "dpi")
~~~
{:.language-r}

Since we do not want to work on all comparisons, we will filter out the samples and conditions that we do not need. Only the mock growth and the _P. syringae_ infected condition will remain.  

~~~
# filter design file (mock versus P. syringae at 7 dpi)
xp_design_mock_vs_infected = xp_design %>% filter(growth == "MgCl2" & dpi == "7")

# Filter count file accordingly (so the rows correspond to the columns of the filtered xp_design file)
counts_filtered = counts[, colnames(counts) %in% xp_design_mock_vs_infected$sample]

## Creation of the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, 
                              colData = xp_design_mock_vs_infected, 
                              design = ~ infected)

~~~
{:.language-r}

It is important to make sure that levels are properly ordered so we are indeed using the _mock_ group as our reference level. A positive gene fold change will for instance signify that the gene is upregulated in the _P. syringae_ condition relatively to the _mock_ condition.  

Please consult [the dedicated section of the DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#factorlevels) on factor levels. 

One way to see how levels are interpreted within the DESeqDataSet object is to display the factor levels. 
~~~
dds$infected
~~~
{:.language-r}

~~~
[1] mock  mock  mock  mock  Pseudomonas_syringae_DC3000
[6] Pseudomonas_syringae_DC3000 Pseudomonas_syringae_DC3000 Pseudomonas_syringae_DC3000
Levels: mock Pseudomonas_syringae_DC3000
~~~
{: .output}

This shows that the _mock_ level comes first before the _Pseudomonas_syringae_DC3000_ level. If this is not correct, you can change it following [the dedicated section of the DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#factorlevels) on factor levels. 


## Differential expression analysis

Differential gene expression analysis will consist of simply two lines of code:
1. The first will call the `DESeq` function on a `DESeqDataSet` object that you've just created under the name `dds`. It will be returned under the same `R` object name `dds`.
2. Then, results are extracted using the `results` function on the `dds` object and results will be extracted as a table under the name `res` (short for results). 

~~~
dds <- DESeq(dds)
res <- results(dds)
~~~
{: .language-r}

The theory beyond DESeq2 differential gene expression analysis is beyond this course but nicely explained [within the DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory). 

> ## Beware of factor levels
> 
> If you do not supply any values to the contrast argument of the `DESeq` function, it will use the first value of the design variable from the design file.
> 
> In our case, we will perform a differential expression analysis between `mock` and `Pseudomonas_syringae_DC3000`. 
> 1. Which of these two is going to be used as the reference level?
> 2. How would you interpret a positive log2 fold change for a given gene?
>
> > ## Solution
> > 1. The `mock` condition is going to be used as the reference level since _m_ from `mock` comes before `P` from `Pseudomonas_syringae_DC3000`.
> > 2. A positive log2 fold change for a gene would mean that this gene is more abundant in `Pseudomonas_syringae_DC3000` than in the `mock` condition.
> {: .solution} 
{: .challenge}

The complete explanation comes from the [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis):
> Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.

A possible preferred way is to specify the comparison of interest explicitly. We are going to name this new result object `res2` and compare it with the previous one called `res`.

~~~
res2 <- results(dds, contrast = c("infected",                       # name of the factor
                                  "Pseudomonas_syringae_DC3000",    # name of the numerator level for fold change
                                  "mock"))                          # name of the denominator level    

~~~
{: .language-r}

If we now compare the `res` and `res2` DESeqResults objects, they should be exactly the same and return a `TRUE` value.

~~~
all_equal(res, res2)
~~~
{: .language-r}










