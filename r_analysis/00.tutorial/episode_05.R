suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))

source("mypca.R")

###############################
# Episode 05: descriptive plots
###############################

## Data import
counts <- read.delim("00.tutorial/counts.txt", header = T, stringsAsFactors = F)
genes <- counts[,1]
counts <- counts[,-1]
row.names(counts) <- genes
xp_design <- read.delim("00.tutorial/experimental_design_modified.txt", header = T, stringsAsFactors = F, colClasses = rep("character",4))

# change col names
colnames(xp_design) <- c("sample", "seed", "infected", "dpi")

# reorder counts columns according to the experimental design file
counts <- counts[,xp_design$sample]

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

# comparison raw and scaled
# only 6 samples taken for clarity
p_raw <- 
  counts[,1:6] %>%
  rownames_to_column("gene") %>% 
  gather(key = "sample", value = "gene_counts", - gene) %>% 
  mutate(gene_counts_log = log10(gene_counts + 1)) %>% 
  ggplot(., aes(x = sample, y = gene_counts_log)) +
  geom_boxplot(fill = "lightgrey") +
  labs(y = "Raw gene counts (log10 transformed)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_scaled <- 
  counts_normalised[,1:6] %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  gather(key = "sample", value = "gene_counts", - gene) %>% 
  mutate(gene_counts_log = log10(gene_counts + 1)) %>% 
  ggplot(., aes(x = sample, y = gene_counts_log)) +
  geom_boxplot(fill = "lightgreen") +
  labs(y = "Scaled gene counts (log10 transformed)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_raw + p_scaled



###################
# PCA: iris dataset
###################

# Loadings refer to parameters like sepal length or petal width
# Scores refer to observations (individual flowers)
dim(iris) # 150 observations, 4 variables



# perform the PCA analysis on only the first 4 variables (skip the Species variable)
iris = t(iris[,1:4]) # transpose the matrix without species character col
pca <- mypca(iris[,1:4], center = TRUE, scale = TRUE)
scores = as.data.frame(pca$scores[,1:2])
scores['Species'] = iris$Species
explained_var = pca$explained_var$exp_var



# plot the scores of the first 2 components
iris_score_plot <- ggplot(scores) +
  geom_point(aes(x = PC1, y = PC2, shape = Species, col = Species), size = 2) +
  xlab(paste0('PC1(',explained_var[1],'%)')) +
  ylab(paste0('PC2(',explained_var[2],'%)')) +
  ggtitle('PCA score plot: PC1 vs PC2')
iris_score_plot
ggsave(filename = "../img/05-iris-score-plot.png", plot = iris_score_plot)


# Exercise
# plot the scores of PC1 and PC3
pca <- mypca(iris[,1:4], center = TRUE, scale = TRUE)
scores = as.data.frame(pca$scores[,1:3])
scores['Species'] = iris$Species
iris_score_plot_pc1_pc3 <- ggplot(scores) +
  geom_point(aes(x = PC1, y = PC3, shape = Species, col = Species), size = 2) +
  xlab(paste0('PC1(',explained_var[1],'%)')) +
  ylab(paste0('PC2(',explained_var[3],'%)')) +
  ggtitle('PCA score plot: PC1 vs PC3')
iris_score_plot_pc1_pc3
ggsave(filename = "../img/05-iris-score-plot-pc3.png", plot = iris_score_plot_pc1_pc3)


# library(reshape2) # to access the melt() function
#
 # reformat the loading data
# loadings <- melt(pca$loadings)
# # rename the columns
# colnames(loadings)<-c("Component","Value")
# # add the 'original' variable names
# loadings['Variable']=as.factor(rep(colnames(iris)[-5],4))
# # plot the loading values per components
# loadings_plot <- ggplot(loadings,
#                         aes(x=Variable,y=Value)) +
#   geom_bar(stat='identity') +
#   facet_wrap(~Component)
# loadings_plot
# ggsave(filename = "../img/05-pca_iris_loadings.png", plot = loadings_plot)
#
# # Plot the first 3 variables with high loading on PC1
# plot_iris_variables <-
#   iris %>%
#   pivot_longer(- Species, names_to = "measurement") %>%
#   ggplot(., aes(x = Species, y = value, fill = Species)) +
#   geom_histogram(stat = "identity") +
#   facet_wrap(~ measurement)
# plot_iris_variables
# ggsave(filename = "../img/05-variables-iris.png", plot = plot_iris_variables)

#####
# VST
#####
# Plot of mean - sd comparison
# Variance - mean plot for all genes
p_mean_sd_scaled <- 
  counts_normalised %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (log10 scale)",
       y = "Gene count standard deviation (log10 scale)") +
  ggtitle("Mean - Standard deviation relationship\n(no variance stabilisation ")
p_mean_sd_scaled
ggsave(filename = "../img/05-mean-sd-gene-before-vst.png", plot = p_mean_sd_scaled)
# while we have corrected for sample to sample differences (e.g. sequencing depth)
# the SD of every gene still has a strong relationship with the mean 
# Our data are "heteroskedastic" meaning that gene variances are not comparable
# we need to address this

# Stabilise the variance to avoid it depending on the mean
dds = estimateDispersions(object = dds, fitType = "parametric", quiet = TRUE)
vsd = varianceStabilizingTransformation(object = dds, 
                                               blind = TRUE, # do not take the formula into account. Best for sample QC
                                               fitType = "parametric")

variance_stabilised_counts <- assay(vsd)

p_mean_sd_vst <- 
  variance_stabilised_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = gene_average, y = gene_stdev)) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (variance stabilised)", 
       y = "Gene count standard deviation (variance stabilised)") +
  ggtitle("Mean - Standard deviation relationship\n(after variance stabilisation ")
p_mean_sd_vst

p_mean_sd_scaled + p_mean_sd_vst
ggsave(filename = "../img/05-comparison-before-after-vst.png")

#################
# PCA: scree plot
#################

pca_results <- mypca(variance_stabilised_counts, 
                     center = TRUE, 
                     scale = TRUE)

# make the plot
ggplot(pca_results$explained_var, 
         aes(x = seq(from = 1, to = nrow(pca_results$explained_var)), 
             y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('Explained variance per component') + 
  geom_bar(stat = "identity") +
  labs(x = "Principal Component number") +
  scale_x_continuous(breaks = seq(
    from = 1, 
    to = nrow(pca_results$explained_var)))
ggsave(filename = "../img/05-screeplot-rnaseq.png")

# 7 components are necessary to catch 50% of the variance
cumsum(pca_results$explained_var) %>% 
  as.data.frame() %>% 
  filter(exp_var > 50) %>% 
  head(n = 1)


#################
# PCA: score plot
#################
scores <- pca_results$scores %>% 
  rownames_to_column("sample")
scores_with_conditions <-  inner_join(x = scores,
                                      y = xp_design,
                                      by = "sample")


# explained variance
# one % variance value per PC
explained_variance <- pca_results$explained_var %>% 
  pull("exp_var")

# non-informative plot
ggplot(scores_with_conditions, 
       aes(PC1, 
           PC2, 
           label = sample)) +
  geom_point(size = 4) +
  geom_text_repel() +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the infection condition")
ggsave(filename = "../img/05-bare-score-plot.png")

# infection score plot
infection_plot <- 
  ggplot(scores_with_conditions, 
       aes(PC1, PC2, color = infected)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the infection condition overlaid")
infection_plot
ggsave(filename = "../img/05-infection-score-plot.png", plot = infection_plot)

# seed score plot
seed_plot <- 
  ggplot(scores_with_conditions, 
       aes(PC1, PC2, color = seed)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the seed condition overlaid")
seed_plot
ggsave(filename = "../img/05-seed-score-plot.png", plot = seed_plot)

# dpi score plot
dpi_plot <- ggplot(scores_with_conditions, 
       aes(PC1, PC2, color = dpi)) +
  geom_point(size = 2) +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the time after infection (dpi) overlaid") +
  scale_color_manual(values = c("orange","blue"))
dpi_plot
ggsave(filename = "../img/05-dpi-score-plot.png", plot = dpi_plot)

# comparing the infection with dpi
infection_plot + dpi_plot
ggsave(filename = "../img/05-infection-dpi-comparison.png", width = 10, height = 5)

# infection + seed score plot
infection_with_seed_plot <- ggplot(scores_with_conditions, 
                   aes(PC1, PC2, color = infected, shape = seed)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the infection and seed inoculation overlaid")
infection_with_seed_plot
ggsave(filename = "../img/05-infection-seed-score-plot.png", plot = infection_with_seed_plot)
