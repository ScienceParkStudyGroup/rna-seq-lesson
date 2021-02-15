suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("DESeq2"))           # Bioconductor package BiocManager::install("DESeq2")
suppressPackageStartupMessages(library("EnhancedVolcano"))  # Bioconductor package BiocManager::install("EnhancedVolcano")
suppressPackageStartupMessages(library("vsn"))
suppressPackageStartupMessages(library("patchwork"))

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

#############
### VST + PCA
#############

# Plot of mean - sd comparison
# Variance - mean plot for all genes
scaled_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") 

# Stabilise the variance to avoid it depending on the mean
dds = estimateDispersions(dds)
vst_counts = getVarianceStabilizedData(object = dds)

#variance_transformed_plot <- meanSdPlot(assay(vst_dds), plot = FALSE, ranks = FALSE)$gg
#not_transformed_plot <- meanSdPlot(assay(dds),ranks = FALSE, plot = FALSE)$gg + scale_x_log10() + scale_y_log10()
#not_transformed_plot
#variance_transformed_plot

# PCA plot on variance transformed dataset
pcaData <- plotPCA(vst(dds), intgroup = c("condition"), returnData = TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot") +
  theme(axis.text = element_text(size = 12))
p
ggsave(filename = "NASA_spaceflight_GLDS38/pca_plot.pdf")


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


##################
# Volcano plot
##################
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


