suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(GGally))

##############
## Data import
##############
counts <- read.delim("counts.txt", header = T, stringsAsFactors = F)
genes <- counts[,1]
counts <- counts[,-1]
row.names(counts) <- genes
xp_design <- read.delim("experimental_design_modified.txt", header = T, stringsAsFactors = F, colClasses = rep("character",4))

# change col names
colnames(xp_design) <- c("sample", "seed", "infected", "dpi")

## Creation of the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = xp_design, 
                              design = ~ seed + infected + dpi)

# normalised counts
dds <- estimateSizeFactors(dds)
counts_normalised = counts(dds, normalized = TRUE)


####################
## sample clustering
####################
# custom colors
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255) # make new colours to extend the 9 colors of the brewer palette (only has 9 color blue levels)



# as we cannot compute each correlation pair 


##################
# diff expression
#################
# Filter design file: keep only mock versus infected
xp_design_mock_vs_infected = xp_design %>% filter(seed == "MgCl2" & dpi == "7") 

# Filter count file accordingly (so the rows correspond to the columns of the filtered xp_design file)
counts_filtered = counts[, colnames(counts) %in% xp_design_mock_vs_infected$sample]

## Creation of the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, 
                              colData = xp_design_mock_vs_infected, 
                              design = ~ infected)



dds <- DESeq(dds)
res <- results(dds, contrast = c("infected",                        # name of the factor
                                  "Pseudomonas_syringae_DC3000",    # name of the numerator level for fold change
                                  "mock"))                          # name of the denominator level 


# how many genes differentially regulated ?
# threshold of p = 0.01
res %>% 
  as.data.frame() %>% 
  filter(padj < 0.01) %>% 
  dim()

# threshold of p = 0.001
res %>% 
  as.data.frame() %>% 
  filter(padj < 0.001) %>% 
  dim()


# one can also see the impact of false discovery rate method
res_wo_fdr = results(dds, 
                     contrast = c("infected", "Pseudomonas_syringae_DC3000", "mock"), 
                     pAdjustMethod = "none")  


#########
# MA plot
#########
# MA plot
plotMA(object = dds, alpha = 0.01)

# shrink effect size
resLFC <- lfcShrink(dds = dds, 
                  res = res,
                  type = "normal",
                  coef = 2) # corresponds to "infected_Pseudomonas_syringae_DC3000_vs_mock" comparison
############## 
# Volcano plot
#############

EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC)
                )

EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 1e-06,
                transcriptPointSize = 2.0,
                FCcutoff = 2, 
                title = "Pseudomonas syringae DC3000 versus mock \n (fold change cutoff = 2, p-value cutoff = 1e-06)",
                legend=c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
                )



######### 
# Heatmap
#########

counts_norm_small <- counts_normalised[1:50,]

# not scaled
pheatmap(counts_norm_small[1:10,], cluster_rows = F, cluster_cols = F)
pheatmap(counts_norm_small[1:20,], cluster_rows = F, cluster_cols = F)
pheatmap(counts_norm_small[1:50,], cluster_rows = F, cluster_cols = F)

# scaled using a log10 transformation
counts_norm_small[counts_norm_small == 0] <- 1

pheatmap(log10(counts_norm_small), cluster_rows = F, cluster_cols = F)


# clearer heatmaps by filtering out genes not differentially expressed
genes_differential = 
  res %>%
  as.data.frame() %>%
  mutate(gene = row.names(res)) %>% 
  filter(padj < 0.01) %>% 
  select(gene) 
  
# use the genes names in res_only_diff to filter the counts_normalised matrix
dim(counts_normalised) # contains all gene info = 33,768 genes

counts_normalised_only_diff = counts_normalised[row.names(counts_normalised) %in% genes_differential$gene, ]

pheatmap(log10(counts_normalised_only_diff + 1), cluster_rows = T, cluster_cols = T)


test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")


annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)), 
  Time = 1:5)
rownames(annotation_col) = paste("Test", 1:10, sep = "")


