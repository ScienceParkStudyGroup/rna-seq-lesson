suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(patchwork))


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
