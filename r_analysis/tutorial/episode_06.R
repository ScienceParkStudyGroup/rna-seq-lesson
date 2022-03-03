suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(apeglm))

source("mypca.R")


#####################################
# Episode 06: differential expression
#####################################

##################################################
# Create the dds object for mock and infected 7dpi
##################################################
# Import libraries
library("DESeq2")
library("tidyverse")

# import the samples to conditions correspodence
xp_design <- read.csv("tutorial/samples_to_conditions.csv", 
                      header = T, 
                      stringsAsFactors = F, 
                      colClasses = rep("character",4))

# filter design file to keep only "mock" and the "infected P. syringae at 7 dpi" conditions.
xp_design_mock_vs_infected = xp_design %>% 
  filter(growth == "MgCl2" & dpi == "7")

# Import the gene raw counts
raw_counts <- read.csv("tutorial/raw_counts.csv", header = T, stringsAsFactors = F) %>% 
  column_to_rownames("Geneid")


# reorder counts columns according to the complete list of samples 
raw_counts <- raw_counts[ , xp_design$sample]

# Filter count file accordingly (to keep only samples present in the filtered xp_design file)
raw_counts_filtered = raw_counts[, colnames(raw_counts) %in% xp_design_mock_vs_infected$sample]

## Creation of the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered, 
                              colData = xp_design_mock_vs_infected, 
                              design = ~ infected)
dds
dds$infected


#########################
# Call differential genes
#########################

dds = DESeq(dds)

##################
# Extract results
#################
all_genes_results = results(dds, 
                            contrast = c("infected",                           # name of the factor
                                         "Pseudomonas_syringae_DC3000",        # name of the numerator level for fold change
                                         "mock"))                               # name of denominator



diff_genes = all_genes_results %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))
head(diff_genes)

################
# Volcano Plot
################

resultsNames(dds)

resLFC <- lfcShrink(dds = dds, 
                    res = all_genes_results,
                    type = "apeglm",
                    coef = "infected_Pseudomonas_syringae_DC3000_vs_mock") # corresponds to "infected_Pseudomonas_syringae_DC3000_vs_mock" comparison