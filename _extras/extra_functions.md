---
title: " Extra functions"
---

# Table of contents

<!-- MarkdownTOC autolink="true" -->

- [1. Custom PCA function](#1-custom-pca-function)
    - [1.1 Credits](#11-credits)
    - [1.2 Usage](#12-usage)
- [2. Custom median of ratios manual normalization](#2-custom-median-of-ratios-manual-normalization)
    - [2.1 Credits](#21-credits)
    - [2.2 Step-by-step procedure](#22-step-by-step-procedure)
        - [Step 1: Upload data](#step-1-upload-data)
        - [Step 2: Take the log of all values](#step-2-take-the-log-of-all-values)
        - [Step 3: Take the average of each row](#step-3-take-the-average-of-each-row)
        - [Step 4: Filter out all of the genes with -Inf as their average](#step-4-filter-out-all-of-the-genes-with--inf-as-their-average)
        - [Step 5: Subtract the gene pseudo-references from log counts](#step-5-subtract-the-gene-pseudo-references-from-log-counts)
        - [Step 6: Find the median of the ratios for each sample](#step-6-find-the-median-of-the-ratios-for-each-sample)
        - [Step 7: Convert medians to scaling factors](#step-7-convert-medians-to-scaling-factors)
        - [Step 8: Divide the original counts \(not log version\) by the scaling factors](#step-8-divide-the-original-counts-not-log-version-by-the-scaling-factors)
    - [2.3 Function to normalize the DESeq2 way](#23-function-to-normalize-the-deseq2-way)
    - [2.4 Section 3: sanity check](#24-section-3-sanity-check)
    - [2.5 References](#25-references)

<!-- /MarkdownTOC -->


# 1. Custom PCA function

## 1.1 Credits
A huge thanks to [Frans van der Kloet](https://bioinformaticslaboratory.eu/team/frans-van-der-kloet/) for his custom function to perform a Principal Component Analysis (PCA). 

This function is used in [episode 05](/05-descriptive-plots/index.html#23-the-iris-data-set).

## 1.2 Usage
Here is the function. To access it in your R environment, two options are available:  
- Option 1: copy and paste the code below in your R console and execute it.   
- Option 2: save the code below as a script called `mypca.R` and source it inside the console `source(mypca.R)` to load the `mypca()` function.  


```bash
# define a custom R function called "mypca()""
mypca <- function(x, center = TRUE, scale = TRUE){
  # Samples should be in rows
  # Variables in the columns

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
```

# 2. Custom median of ratios manual normalization

## 2.1 Credits 
A huge thanks to [Jolanta Szkodon](https://github.com/jszkodon) for her custom function to reproduce DESeq2 normalization. 

The median of ratios method is a normalization method used to account for sequencing differences in library size and RNA composition of samples. This makes it a suitable method to use when comparing between samples. It does not account for differences in gene length, making it unsuitable for within sample comparisons.

The median of ratios method is applied in the DESeq2 package, used for differential expression testing. In the DESeq2 package, normalization is conducted easily with a function. Although this is convenient, it is useful to understand the step-by-step breakdown of the median of ratios method. 

Sample expression count data was taken from the [GSE49110 dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49110) and cut to create a manageable dataset for the tutorial. 

In [Section 1](#section-1-step-by-step-procedure), steps of the median of ratios method are shown.  
In [Section 2](#section-two-a-function-to-normalize-the-deseq2-way), a function including all steps is shown.   
In [Section 3](#section-3-sanity-check), the normalized counts are compared with the DESeq2 method to ensure the results are the same. References used are listed at the end. 

## 2.2 Step-by-step procedure 

### Step 1: Upload data 
Rows (genes) and columns (samples) are seen below. 
~~~
data = read.csv("sample_count_data.csv", row.names = 1)
head(data)
~~~
{: .language-r}

~~~
          siC1 siC2 siE11 siE12 siE21 siE22 siE31 siE32
12SBSRNA4    0    2     3     3     2     4     6     2
A1BG        50   52    86    80    66    94    67    83
A1BG-AS1     7    9    12    10    11    18     7     8
A1CF         2    1     4     7     7     1     1     7
A2LD1       58   76    57    52    67    97   101    98
A2M          0    1     0     1     1     0     0     0
~~~
{: .output}

The following few steps are used to create a pseudo-reference sample for each gene. To do this, we need to calculate the geometric mean for each gene. 

### Step 2: Take the log of all values
~~~
log_data = log(data)
head(log_data)
~~~
{: .language-r}

~~~
               siC1      siC2    siE11    siE12     siE21    siE22    siE31     siE32
12SBSRNA4      -Inf 0.6931472 1.098612 1.098612 0.6931472 1.386294 1.791759 0.6931472
A1BG      3.9120230 3.9512437 4.454347 4.382027 4.1896547 4.543295 4.204693 4.4188406
A1BG-AS1  1.9459101 2.1972246 2.484907 2.302585 2.3978953 2.890372 1.945910 2.0794415
A1CF      0.6931472 0.0000000 1.386294 1.945910 1.9459101 0.000000 0.000000 1.9459101
A2LD1     4.0604430 4.3307333 4.043051 3.951244 4.2046926 4.574711 4.615121 4.5849675
A2M            -Inf 0.0000000     -Inf 0.000000 0.0000000     -Inf     -Inf      -Inf
~~~
{: .output}

### Step 3: Take the average of each row
~~~
library(dplyr)
library(tibble)

log_data = log_data %>% 
             rownames_to_column('gene') %>% 
             mutate (pseudo_reference = rowMeans(log_data))

head(log_data)
~~~
{: .language-r}

~~~
       gene      siC1      siC2    siE11    siE12     siE21    siE22    siE31     siE32 pseudo_reference
1 12SBSRNA4      -Inf 0.6931472 1.098612 1.098612 0.6931472 1.386294 1.791759 0.6931472             -Inf
2      A1BG 3.9120230 3.9512437 4.454347 4.382027 4.1896547 4.543295 4.204693 4.4188406        4.2570154
3  A1BG-AS1 1.9459101 2.1972246 2.484907 2.302585 2.3978953 2.890372 1.945910 2.0794415        2.2805306
4      A1CF 0.6931472 0.0000000 1.386294 1.945910 1.9459101 0.000000 0.000000 1.9459101        0.9896465
5     A2LD1 4.0604430 4.3307333 4.043051 3.951244 4.2046926 4.574711 4.615121 4.5849675        4.2956204
6       A2M      -Inf 0.0000000     -Inf 0.000000 0.0000000     -Inf     -Inf      -Inf             -Inf

~~~
{: .output}


### Step 4: Filter out all of the genes with -Inf as their average
They will not be used to calculate the median.
~~~
filtered_log_data = log_data %>% filter(pseudo_reference != "-Inf")
head(filtered_log_data)
~~~
{: .language-r}

~~~
      gene      siC1     siC2    siE11    siE12    siE21    siE22     siE31    siE32 pseudo_reference
1     A1BG 3.9120230 3.951244 4.454347 4.382027 4.189655 4.543295 4.2046926 4.418841        4.2570154
2 A1BG-AS1 1.9459101 2.197225 2.484907 2.302585 2.397895 2.890372 1.9459101 2.079442        2.2805306
3     A1CF 0.6931472 0.000000 1.386294 1.945910 1.945910 0.000000 0.0000000 1.945910        0.9896465
4    A2LD1 4.0604430 4.330733 4.043051 3.951244 4.204693 4.574711 4.6151205 4.584967        4.2956204
5   A4GALT 6.1506028 6.100319 6.270988 6.565265 6.335054 6.428105 6.3952616 6.171701        6.3021621
6    A4GNT 0.6931472 0.000000 0.000000 1.609438 1.945910 1.098612 0.6931472 0.000000        0.7550318
~~~
{: .output}

### Step 5: Subtract the gene pseudo-references from log counts
In this step, you are subtracting the average of the logs from the log of the original data. Log(x/y) is a ratio. 
~~~
ratio_data = sweep(filtered_log_data[,2:9], 1, filtered_log_data[,10], "-")
head(ratio_data)
~~~
{: .language-r}

~~~
         siC1        siC2       siE11       siE12       siE21      siE22       siE31      siE32
1 -0.34499242 -0.30577171  0.19733187  0.12501121 -0.06736068  0.2862794 -0.05232281  0.1618252
2 -0.33462050 -0.08330607  0.20437600  0.02205444  0.11736462  0.6098411 -0.33462050 -0.2010891
3 -0.29649932 -0.98964650  0.39664786  0.95626365  0.95626365 -0.9896465 -0.98964650  0.9562637
4 -0.23517736  0.03511297 -0.25256910 -0.34437665 -0.09092775  0.2790906  0.31950015  0.2893471
5 -0.15155934 -0.20184315 -0.03117367  0.26310286  0.03289215  0.1259432  0.09309949 -0.1304615
6 -0.06188466 -0.75503184 -0.75503184  0.85440607  1.19087831  0.3435804 -0.06188466 -0.7550318
~~~
{: .output}

### Step 6: Find the median of the ratios for each sample
~~~
sample_medians = apply(ratio_data, 2, median)
sample_medians
~~~
{: .language-r}

~~~
       siC1        siC2       siE11       siE12       siE21       siE22       siE31       siE32 
-0.02605145 -0.04549215 -0.07749397  0.12030642  0.06268842  0.05039392 -0.06025587  0.03220589 
~~~
{: .output}

### Step 7: Convert medians to scaling factors
~~~
# scaling factor = e^median
scaling_factors = exp(sample_medians)
scaling_factors
~~~
{: .language-r}

~~~
     siC1      siC2     siE11     siE12     siE21     siE22     siE31     siE32 
0.9742850 0.9555271 0.9254326 1.1278424 1.0646951 1.0516853 0.9415236 1.0327301 
~~~
{: .output}

### Step 8: Divide the original counts (not log version) by the scaling factors 
~~~
manually_normalized = sweep(data, 2, scaling_factors, "/")
head(manually_normalized)
~~~
{: .language-r}

~~~
               siC1      siC2     siE11      siE12      siE21      siE22      siE31     siE32
12SBSRNA4  0.000000  2.093086  3.241727  2.6599461  1.8784721  3.8034192   6.372650  1.936614
A1BG      51.319688 54.420225 92.929511 70.9318970 61.9895807 89.3803503  71.161255 80.369498
A1BG-AS1   7.184756  9.418885 12.966909  8.8664871 10.3315968 17.1153862   7.434758  7.746458
A1CF       2.052788  1.046543  4.322303  6.2065410  6.5746525  0.9508548   1.062108  6.778150
A2LD1     59.530838 79.537252 61.592816 46.1057330 62.9288167 92.2329146 107.272936 94.894106
A2M        0.000000  1.046543  0.000000  0.8866487  0.9392361  0.0000000   0.000000  0.000000
~~~
{: .output}

## 2.3 Function to normalize the DESeq2 way

Below is one function to "manually" normalize data with the median of ratios method. Import data before using the function.

**mor** stands for "median of ratios".

~~~
mor_normalization = function(data){
  library(dplyr)
  library(tibble)

  # take the log
  log_data = log(data) 
  
  # find the psuedo-references per sample by taking the geometric mean
  log_data = log_data %>% 
               rownames_to_column('gene') %>% 
               mutate (gene_averages = rowMeans(log_data)) %>% 
               filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  return(manually_normalized)
}
~~~
{: .language-r}


You can check that this function returns the same result as the step-by-step procedure above :arrow_up: . 

~~~
mor_normalization(data) == manually_normalized
~~~
{: .language-r}

~~~
          siC1 siC2 siE11 siE12 siE21 siE22 siE31 siE32
12SBSRNA4 TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A1BG      TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A1BG-AS1  TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A1CF      TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A2LD1     TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
... (more lines) ...
~~~
{: .output}

## 2.4 Section 3: sanity check

A sanity check to see if the normalization within DESeq2 is identical to our manually normalized data. 

~~~
library(DESeq2)
~~~
{: .language-r}

~~~
# samples (columns names) of the data should be named
samples = as.data.frame(colnames(data))

# create a DESeqDataSet object. The design can be altered based on experimental design. A design of 1 means no design. 
dds = DESeqDataSetFromMatrix(countData = data, colData = samples, design = ~1)

# this function generates the size factors
dds = estimateSizeFactors(dds)
~~~
{: .language-r}

We can compare the scaling factors that were found manually vs. the scaling factors (termed, "size factors") from DESeq2. 
~~~
# scaling_factors were manually computed using our mor_normalization function
# sizeFactors(dds) is used to find the scaling factors from DESeq2
scaling_factors == sizeFactors(dds)
~~~
{: .language-r}

~~~
 siC1  siC2 siE11 siE12 siE21 siE22 siE31 siE32 
 TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
~~~
{: .output}

Lastly, we can use `counts()` to normalize the original data in `DESeq2`. When `counts()` is used with the `normalized = TRUE` arguments, the counts will be normalized. When you do not use this argument, counts(dds) returns the original counts. 
~~~
normalized_deseq2 = counts(dds, normalized = TRUE)
~~~
{: .language-r}

The normalized counts are identical between DESeq2 and the manual method. 
~~~
normalized_deseq2 == manually_normalized
~~~
{: .language-r}

~~~
          siC1 siC2 siE11 siE12 siE21 siE22 siE31 siE32
12SBSRNA4 TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A1BG      TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A1BG-AS1  TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A1CF      TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A2LD1     TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A2M       TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A2ML1     TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
A2MP1     TRUE TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
... (more lines) ...
~~~
{: .output}

## 2.5 References
1. Love et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. [doi:10.1186/s13059-014-0550-8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).
2. [Statquest](https://www.youtube.com/watch?v=UFB993xufUU)
3. [DGE count normalization](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)



