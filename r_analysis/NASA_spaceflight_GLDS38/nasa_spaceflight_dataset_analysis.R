suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("DESeq2"))           # Bioconductor package BiocManager::install("DESeq2")
suppressPackageStartupMessages(library("EnhancedVolcano"))  # Bioconductor package BiocManager::install("EnhancedVolcano")
suppressPackageStartupMessages(library("vsn"))
suppressPackageStartupMessages(library("patchwork"))

source("mypca.R")

new_sample_names = c(
  paste0("FLT", seq(1,3,1)),
  paste0("GC", seq(1,3,1))
)

###############
# External info
###############
arabidopsis2uniprot <- read_tsv(file = "NASA_spaceflight_GLDS38/Arabidopsis_to_Uniprot_idmapping.txt", 
                    col_names = c("uniprot_id","gene")) %>% 
  group_by(gene) %>% 
  slice_head()   %>%                                # to take only one value if duplicated Uniprot ids
  mutate(ipath_id = paste0("UNIPROT:", uniprot_id)) # for compatibility with iPath https://pathways.embl.de/

   

#######################
# Create DESeq2 dataset
#######################
raw_counts <- read.csv("NASA_spaceflight_GLDS38/GLDS38_raw_counts.csv", 
                       stringsAsFactors = F) %>% 
  column_to_rownames("gene")

sample2condition <- read.csv("NASA_spaceflight_GLDS38/samples_to_condition.csv",
                             stringsAsFactors = F)


dds <- DESeqDataSetFromMatrix(countData = raw_counts, 
                              colData = sample2condition, 
                              design = ~ condition)


###################
# Get scaled counts
###################
dds <- estimateSizeFactors(dds)
scaled_counts = counts(dds, normalized = TRUE) 
write.csv(x = scaled_counts, 
          file = "NASA_spaceflight_GLDS38/GLDS38_scaled_counts.csv",
          quote = F,
          row.names = TRUE)


as.data.frame(sizeFactors(dds)) %>% 
  rownames_to_column("sample") %>% 
  dplyr::rename("size_factor" = "sizeFactors(dds)") %>% 
  ggplot(., aes(x = sample, y = size_factor)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Calculated size factor (correction factor)") +
  geom_hline(yintercept = 1, colour = "red")

ggsave(filename = "NASA_spaceflight_GLDS38/size_factor_plot.pdf", width = 10, height = 7) 

###############################
# Compare raw and scaled counts
###############################
colnames(raw_counts) = new_sample_names
colnames(scaled_counts) = new_sample_names

sample2condition$sample = new_sample_names

p_raw <- 
  raw_counts %>%
  rownames_to_column("gene") %>% 
  gather(key = "sample", value = "gene_counts", - gene) %>% 
  mutate(gene_counts_log = log10(gene_counts + 1)) %>% 
  ggplot(., aes(x = sample, y = gene_counts_log)) +
  geom_boxplot(fill = "lightgrey") +
  labs(y = "Raw gene counts (log10 transformed)")
  
p_scaled <- 
  scaled_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  gather(key = "sample", value = "gene_counts", - gene) %>% 
  mutate(gene_counts_log = log10(gene_counts + 1)) %>% 
  ggplot(., aes(x = sample, y = gene_counts_log)) +
  geom_boxplot(fill = "lightgreen") +
  labs(y = "Scaled gene counts (log10 transformed)")

p_raw + p_scaled


#############
### VST + PCA
#############
# Stabilise the variance to avoid it depending on the mean
dds = estimateDispersions(object = dds, fitType = "mean")
vst_counts = getVarianceStabilizedData(object = dds)

colnames(vst_counts) = new_sample_names

# Plot of mean - sd comparison
# Variance - mean plot for all genes
p_mean_sd_scaled <- 
  scaled_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") 
p_mean_sd_scaled

p_mean_sd_vst <- 
  vst_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = gene_average, y = gene_stdev)) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") 
p_mean_sd_vst

p_mean_sd_scaled + p_mean_sd_vst


# PCA plot with the mypca() on scaled counts
scaled_counts <- t(scaled_counts) # samples in rows and genes in columns
pca_results <- mypca(scaled_counts, scale = TRUE)
percentage_variance <- as.data.frame(pca_results$explained_var)

# make the plot
scree_plot <- 
  percentage_variance %>% 
  rownames_to_column("PC") %>% 
  ggplot(., aes(x = PC, y = exp_var)) +
    ylab('explained variance (%)') + 
    ggtitle('explained variance per component') + 
    geom_bar(stat = "identity")
scree_plot

scores <- pca_results$scores %>% 
  rownames_to_column("sample") %>% 
  left_join(., y = sample2condition, by = "sample")

p_pca_scaled <- ggplot(scores, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentage_variance[1,],"% variance")) +
  ylab(paste0("PC2: ",percentage_variance[2,],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot (scaled counts") +
  theme(axis.text = element_text(size = 12))
p_pca_scaled
ggsave(filename = "NASA_spaceflight_GLDS38/pca_plot_scaled_counts.pdf")


# PCA on the vst counts
t_vst_counts <- t(vst_counts)
pca_results_vst <- mypca(t_vst_counts)
scores_vst <- pca_results_vst$scores %>% 
  rownames_to_column("sample") %>% 
  left_join(., y = sample2condition, by = "sample")
percentage_variance_vst <- as.data.frame(pca_results_vst$explained_var)

p_vst_counts <- ggplot(scores_vst, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentage_variance_vst[1,],"% variance")) +
  ylab(paste0("PC2: ",percentage_variance_vst[2,],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot (vst transformed") +
  theme(axis.text = element_text(size = 12))
p_vst_counts

p_pca_scaled + p_vst_counts

###################
### diff expression
###################

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "space_flight", "ground_control")) %>% 
  as.data.frame()

# histogram of raw and adjusted p-values
plot_raw_pvalues <- ggplot(res, aes(x = pvalue)) +
  geom_histogram(fill = "grey") + ggtitle("Raw p-values") + labs(y = "Number of occurences")
plot_adj_pvalues <- ggplot(res, aes(x = padj)) +
  geom_histogram(fill = "lightblue") + ggtitle("Adjusted p-values") + labs(y = "Number of occurences")
plot_raw_pvalues + plot_adj_pvalues
ggsave(filename = "NASA_spaceflight_GLDS38/pvalues_histograms.pdf")

# distribution of log2 fold changes
res %>% 
  drop_na() %>% 
  ggplot(., aes(x = log2FoldChange)) +
  geom_histogram() +
  geom_vline(xintercept = 1, colour= "red") +
  geom_vline(xintercept = -1, colour= "blue") +
  scale_x_continuous(breaks = seq(-5,+5,1), limits = c(-5,+5)) +
  theme(text = element_text(size = 20))
ggsave(filename = "NASA_spaceflight_GLDS38/distri_log2fc.pdf")


calculate_percentile_of_log2fc_value <- res %>% 
  drop_na() %>% 
  with(., ecdf(log2FoldChange)) # empirical cumulative distribution function 

# calculate the nth percentile of selected log2 FC
calculate_percentile_of_log2fc_value(-2) * 100 %>%  round(digits = 2) # -> the 6th percentile 
calculate_percentile_of_log2fc_value(+2)* 100 %>%  round(digits = 2)  # -> the 96th percentile  

# genes upregulated by microgravity
res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%  
  as_tibble() %>% 
  filter(padj < 0.01) %>% 
  filter(log2FoldChange > 1) %>% 
  left_join(x = ., y = arabidopsis2uniprot, by = "gene") %>% 
  write_csv(file = "NASA_spaceflight_GLDS38/pos_diff_genes.csv")

# genes downregulated by microgravity
res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%  
  as_tibble() %>% 
  filter(padj < 0.01) %>% 
  filter(log2FoldChange < -1) %>% 
  left_join(x = ., y = arabidopsis2uniprot, by = "gene") %>% 
  write_csv(file = "NASA_spaceflight_GLDS38/neg_diff_genes.csv")

# Top 20 genes (absolute log2FC)
# genes upregulated by microgravity
res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%  
  as_tibble() %>% 
  filter(padj < 0.01) %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>% 
  arrange(desc(abs_log2FoldChange)) %>% 
  top_n(n = 20) %>% 
  filter(log2FoldChange > 0)

res %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>%  
  as_tibble() %>% 
  filter(gene == "AT2G33380") 

################# Volcano plot #################
res <- results(dds, contrast = c("condition", "space_flight", "ground_control"))
resLFC <- lfcShrink(dds = dds, 
                    res = res,
                    type = "normal",
                    coef = 2) # corresponds to "infected_Pseudomonas_syringae_DC3000_vs_mock" comparison



EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-5, +5),
                ylim = c(0,50),
                pCutoff = 0.01,
                transcriptPointSize = 2.0,
                FCcutoff = 1, 
                title = "Volcano plot",
                legend=c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')) + 
  guides(legend = NULL)


#########################################################
# Extra: mypca(scale = TRUE) versus mypca(scale = FALSE)
########################################################

# Scaled counts (NOT variance stabilised)
pca1 <- mypca(scaled_counts, scale = TRUE)
pca2 <- mypca(scaled_counts, scale = FALSE)
pca3 <- mypca(t_vst_counts, scale = TRUE)
pca4 <- mypca(scaled_counts, scale = FALSE)

scores1 <- pca1$scores %>% 
  rownames_to_column("sample") %>% 
  left_join(., y = sample2condition, by = "sample")
percentage_variance1 <- as.data.frame(pca1$explained_var)

p1 <- ggplot(scores1, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentage_variance1[1,],"% variance")) +
  ylab(paste0("PC2: ",percentage_variance1[2,],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot on scaled counts + pca scaling TRUE") +
  theme(axis.text = element_text(size = 12))
p1

scores2 <- pca2$scores %>% 
  rownames_to_column("sample") %>% 
  left_join(., y = sample2condition, by = "sample")
percentage_variance2 <- as.data.frame(pca2$explained_var)
p2 <- ggplot(scores2, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentage_variance2[1,],"% variance")) +
  ylab(paste0("PC2: ",percentage_variance2[2,],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot (vst transformed") +
  ggtitle("PCA plot on scaled counts + pca scaling FALSE") +
  theme(axis.text = element_text(size = 12))
p2

scores3 <- pca3$scores %>% 
  rownames_to_column("sample") %>% 
  left_join(., y = sample2condition, by = "sample")
percentage_variance3 <- as.data.frame(pca3$explained_var)
p3 <- ggplot(scores3, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentage_variance3[1,],"% variance")) +
  ylab(paste0("PC2: ",percentage_variance3[2,],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot (vst transformed") +
  ggtitle("PCA plot on variance stabilised counts + pca scaling TRUE") +
  theme(axis.text = element_text(size = 12))
p3

scores4 <- pca4$scores %>% 
  rownames_to_column("sample") %>% 
  left_join(., y = sample2condition, by = "sample")
percentage_variance4 <- as.data.frame(pca4$explained_var)
p4 <- ggplot(scores4, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentage_variance4[1,],"% variance")) +
  ylab(paste0("PC2: ",percentage_variance4[2,],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot (vst transformed") +
  ggtitle("PCA plot on variance stabilised counts + pca scaling FALSE") +
  theme(axis.text = element_text(size = 12))
p4

p1 + p2 + p3 + p4

bind_cols(percentage_variance1,
          percentage_variance2,
          percentage_variance3, 
          percentage_variance4) %>% 
  rownames_to_column("PC")

ggsave(filename = "NASA_spaceflight_GLDS38/plots/effects_of_vst_versus_pca_scaling.png", 
       width = 12, 
       height = 8)
