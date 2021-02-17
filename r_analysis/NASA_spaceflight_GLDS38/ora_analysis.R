library("biomartr")
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("org.At.tair.db"))
library("biomaRt")  # only use to remove cache bug

diff_genes <- read.csv("NASA_spaceflight_GLDS38/pos_diff_genes.csv", stringsAsFactors = F)

all_arabidopsis_genes <- read.csv("NASA_spaceflight_GLDS38/GLDS38_raw_counts.csv", 
                                    header = T, 
                                    stringsAsFactors = F)[,1] # directly selects the gene column

attributes_to_retrieve = c("tair_symbol", "uniprotswissprot","entrezgene_id") 

all_arabidopsis_genes_annotated <- biomartr::biomart(genes = all_arabidopsis_genes,
                                                     mart       = "plants_mart",                 
                                                     dataset    = "athaliana_eg_gene",           
                                                     attributes = attributes_to_retrieve,        
                                                     filters =  "ensembl_gene_id" )  
all_arabidopsis_genes_annotated$entrezgene_id = as.character(
  all_arabidopsis_genes_annotated$entrezgene_id) 

diff_arabidopsis_genes_annotated <- biomartr::biomart(genes = diff_genes$gene,
                                                      mart       = "plants_mart",                 
                                                      dataset    = "athaliana_eg_gene",           
                                                      attributes = attributes_to_retrieve,        
                                                      filters =  "ensembl_gene_id" )  


###################
# ORA gene ontology
###################

# BP
ora_analysis_bp <- enrichGO(gene = diff_arabidopsis_genes_annotated$entrezgene_id, 
                            universe = all_arabidopsis_genes_annotated$entrezgene_id, 
                            OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "ENTREZID",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE, 
                            pool = FALSE) 
ora_analysis_simplified_bp <- simplify(ora_analysis_bp)

# CC
ora_analysis_cc <- enrichGO(gene = diff_arabidopsis_genes_annotated$entrezgene_id, 
                            universe = all_arabidopsis_genes_annotated$entrezgene_id, 
                            OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "ENTREZID",
                            ont = "CC",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE, 
                            pool = FALSE) 
ora_analysis_simplified_cc <- simplify(ora_analysis_cc)

# MF
ora_analysis_mf <- enrichGO(gene = diff_arabidopsis_genes_annotated$entrezgene_id, 
                            universe = all_arabidopsis_genes_annotated$entrezgene_id, 
                            OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "ENTREZID",
                            ont = "MF",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE, 
                            pool = FALSE) 

ora_analysis_simplified_mf <- simplify(ora_analysis_mf)

rm(ora_analysis_bp)
rm(ora_analysis_mf)
rm(ora_analysis_cc)
#####################################################
# Gene Set Enrichment Analysis Gene Ontology analysis
#####################################################

# Preparation of geneList object
# Column 1: gene ID
# Column 2: numeric value (log2 fold change)
# Column 3: 

geneList = 
  diff_genes %>% 
  dplyr::select(gene, log2FoldChange) %>% 
  dplyr::arrange(desc(log2FoldChange)) %>% 
  with(., set_names(log2FoldChange, gene))
  

gsea_analysis_bp <- gseGO(geneList = geneList,
                          OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                          keyType = "ENTREZID",
                          ont = "BP",
                          nPerm = 1000,
                          minGSSize = 5,
                          maxGSSize = 10*length(geneList), # 10% of input
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          by = 'fgsea')
                          
