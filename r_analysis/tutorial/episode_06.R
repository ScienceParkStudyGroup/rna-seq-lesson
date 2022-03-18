suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(apeglm))

source("tutorial/custom_functions/mypca.R")
source("tutorial/custom_functions/mor_normalization.R")


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
                    coef = "infected_Pseudomonas_syringae_DC3000_vs_mock") 


#########
# Heatmap
#########

#############
# Scale counts
##############
scaled_counts <- mor_normalization(raw_counts_filtered)
head(scaled_counts)

# Compare scaled and raw
long_scaled_counts = 
  scaled_counts %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>% 
  mutate(scaled = "yes")

long_raw_counts = 
  raw_counts %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>% 
  mutate(scaled = "no")

long_raw_and_scaled_counts = bind_rows(long_raw_counts, long_scaled_counts)

ggplot(long_raw_and_scaled_counts, 
       aes(x = scaled, y = counts + 1, fill = scaled)) +
  geom_violin() +
  scale_y_log10() +
  labs(x = "Gene counts scaled/normalised?", y = "Gene counts (raw or scaled)")

##############################
# Keep only differential genes
##############################

mat_for_heatmap <- 
  raw_counts_filtered %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% diff_genes$genes) %>%  # diff_genes contains genes with padj < 0.01 
  column_to_rownames("gene")

counts_normalised_only_diff_genes = 
  mor_normalization(raw_counts) %>%             # normalize the counts using our custom function
  rownames_to_column("genes") %>%               
  pivot_longer(- genes,                         
               names_to = "sample", 
               values_to = "counts") %>% 
  filter(genes %in% diff_genes$genes) %>%             
  pivot_wider(names_from = "sample",            
              values_from = "counts")  %>%      
  column_to_rownames("genes")                   # the gene column is converted back to row names to create a matrix usable with pheatmap

dim(counts_normalised_only_diff_genes) 

pheatmap(counts_normalised_only_diff_genes, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "none",
         show_rownames = FALSE, 
         show_colnames = TRUE)
counts_scaled = 
  counts_normalised_only_diff_genes %>% 
  t(.) %>%                              # transpose to have the genes in columns 
  scale() %>%                           # scale(x, center = TRUE, scale = TRUE) 
  t(.)                                  # back in original shape

pheatmap(counts_scaled, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         scale = "none",            # already done "manually"
         show_colnames = TRUE)

pheatmap(counts_scaled, 
         cluster_rows = TRUE,                      
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         main = "Clustering on")

counts_scaled_filtered = 
  counts_scaled %>% 
  as.data.frame() %>%
  dplyr::select(xp_design_mock_vs_infected$sample) # keep the 8 samples

anno_col_info = xp_design_mock_vs_infected %>% column_to_rownames("sample")

anno_info_colors = list(
  seed = c(MgCl2 = "#d8b365"),
  infected = c(mock = "lightgrey", 
               Pseudomonas_syringae_DC3000 = "black"),
  dpi = c("7" = "dodgerblue4")
)

pheatmap(counts_scaled_filtered, 
         cluster_rows = TRUE,                       
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         annotation_col = anno_col_info,
         annotation_colors = anno_info_colors,
         main = "Clustering with ward method")
