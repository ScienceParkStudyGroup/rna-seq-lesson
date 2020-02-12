suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))


###############################
# Episode 05: descriptive plots
###############################

## Data import
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


#################################
## Size factors and normalisation
#################################
# normalised counts
dds <- estimateSizeFactors(dds)

size_factors_df <- data.frame(sample = names(sizeFactors(dds)), 
                              size = sizeFactors(dds))

# add the experimental condition of interest for plot labelling
size_factors_df <- left_join(size_factors_df, xp_design, by = "sample")

# sort by seed condition and by infected condition
size_factors_df <- size_factors_df %>% 
  arrange(seed, infected)

size_factors_df$sample = factor(size_factors_df$sample, levels = size_factors_df$sample)

# plot
ggplot(size_factors_df, aes(x = sample, y = size, colour = seed)) + 
  geom_segment(aes(x = sample, xend = sample, y = 0, yend = size), color="grey") +
  geom_point(size = 4) + 
  coord_flip() +
  theme_grey() +
  facet_wrap(~ infected) +
  scale_y_continuous(limits = c(0,2))

# extract normalised counts
counts_normalised = counts(dds, normalized = TRUE)

# we reorder samples according to their experimental conditions
counts_normalised = counts[, xp_design$sample]

####################
## sample correlations
####################

# self correlation
pairs(x = counts_normalised[,c(1,1)],pch = 19, log = "xy")
cor(counts_normalised[,c(1,1)])

# correlation between biological replicates
pairs(x = counts_normalised[,1:4],pch = 19, log = "xy")
cor(counts_normalised[,c(1:4)])

# weakly correlated samples 
pairs(x = counts_normalised[,c(1,9,17,25)],pch = 19, log = "xy")
cor(counts_normalised[,c(1,9,17,25)])

#counts_normalised_scaled = t(scale(t(counts_normalised), center = T, scale = F))
#ggcorrplot::ggcorrplot(cor(counts_normalised_scaled))


#counts_normalised_scaled = t(scale(t(log10(counts_normalised+1)), center = T, scale = F))
#ggcorrplot::ggcorrplot(cor(counts_normalised_scaled))

##########
# PCA plot
##########
# The PCA plot, log transform the data to prevent that the largest values will dominate the analyses. Center the data in order to be to analyze the relation between the objects and not the relation of the objects between the origin (0).

counts_norm_trans = t(scale(t(log10(counts_normalised + 1)),scale = FALSE, center = TRUE))

# perform the PCA analysis
pca <- princomp(counts_norm_trans)

# screeplot
screeplot(pca, ylim=c(0,0.25))

# copy the data from the first 2 components
T12 = as.data.frame(pca$loadings[,1:2])
# add the experimental design data to the PCA data.frame
T12c <- cbind.data.frame(T12,xp_design)
# calculate the explained variance per component
explained_var = round(pca$sdev^2/sum(pca$sdev^2)*100,1)
# plot the PCA score plot

p <- ggplot(data = T12c, 
            aes(x = Comp.1, y = Comp.2, col = seed, shape = infected, size = dpi)
            ) + 
  geom_point() + 
  xlab(paste0('PC1(',explained_var[1],'%)')) + 
  ylab(paste0('PC2(',explained_var[2],'%)')) + 
  ggtitle('PCA log10 transformed (centered) data')


# PCA using the plotPCA function
# variance-stabilizing transformation
vst_dds <- vst(dds)

# plot the PCA
# without customisation
plotPCA(vst_dds, intgroup = c("seed", "infected", "dpi"))

# customised PCA plot
pcaData <- plotPCA(vst_dds, intgroup = c("seed", "infected", "dpi"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = seed, shape = infected, size = dpi)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA customized plot")


#####################################
# Episode 06: differential expression
#####################################

##################
# diff expression
#################
# Filter design file: keep only mock versus infected
xp_design_mock_vs_infected = xp_design %>% filter(seed == "MgCl2" & dpi == "7") 

# Filter count file accordingly (so the rows correspond to the columns of the filtered xp_design file)
counts_filtered = counts[, colnames(counts) %in% xp_design_mock_vs_infected$sample]

## Creation of the DESeqDataSet object
dds2 <- DESeqDataSetFromMatrix(countData = counts_filtered, 
                              colData = xp_design_mock_vs_infected, 
                              design = ~ infected)



dds2 <- DESeq(dds2)
res <- results(dds2, contrast = c("infected",                        # name of the factor
                                  "Pseudomonas_syringae_DC3000",     # name of the numerator level for fold change
                                  "mock"))                           # name of the denominator level 


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
res_wo_fdr = results(dds2, 
                     contrast = c("infected", "Pseudomonas_syringae_DC3000", "mock"), 
                     pAdjustMethod = "none")  


#########
# MA plot
#########
# MA plot
plotMA(object = dds2, alpha = 0.01)

# shrink effect size
resLFC <- lfcShrink(dds = dds2, 
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
counts_normalised_scaled = t(scale(t(counts_normalised), center = T, scale = F))

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


