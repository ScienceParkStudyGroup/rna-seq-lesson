---
title: " Custom PCA function"
---

A huge thanks to [Frans van der Kloet](https://bioinformaticslaboratory.eu/team/frans-van-der-kloet/) for her custom function to perform Principal Component Analysis. 

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
