if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("ExpressionAtlas")
suppressMessages(library("ExpressionAtlas"))
suppressMessages(library("tidyverse"))

# search and retrive ArrayExpress accessions relatedo to keywords
atlasRes <- searchAtlasExperiments( properties = "salt", species = "rice") 
rnaseqExps <- getAtlasData(
  experimentAccessions = atlasRes$Accession[grep("rna-seq",atlasRes$Type,ignore.case = TRUE)] 
)

################
# one experiment
################
mtab1625 <- allExps[[ "E-MTAB-1625" ]]

# colData
colData(mtab1625$rnaseq)
# raw counts
cts <- assays(mtab1625$rnaseq)$counts
cts[1:5,1:5]



##########
# Tomato
#########

# search and retrive ArrayExpress accessions related to keywords
atlasRes <- searchAtlasExperiments( properties = "tissue", species = "tomato") 
atlasRes
rnaseqExps <- getAtlasData(
  experimentAccessions = atlasRes$Accession[grep("rna-seq",atlasRes$Type,ignore.case = TRUE)] 
)
mtab4812 <- rnaseqExps[["E-MTAB-4812"]]
colData(mtab4812$rnaseq)
cts <- assays(mtab4812$rnaseq)$counts
cts[1:5,1:5]
