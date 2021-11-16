suppressPackageStartupMessages(library("biomartr"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("org.At.tair.db"))
suppressPackageStartupMessages(library("biomaRt"))  # only use to remove cache bug
suppressPackageStartupMessages(library("splitstackshape"))

####################
# 02: Annotate genes
####################
diff_genes <- read_delim("00.tutorial/differential_genes.tsv", delim = "\t")

biomartr::organismBM(organism = "Arabidopsis thaliana")

arabido_attributes = 
  biomartr::organismAttributes("Arabidopsis thaliana") %>% 
  filter(dataset == "athaliana_eg_gene")
arabido_attributes

attributes_to_retrieve = c("tair_symbol", "entrezgene_id")

result_BM <- biomartr::biomart( genes      = diff_genes$genes,                  # genes were retrieved using biomartr::getGenome()
                                mart       = "plants_mart",                     # marts were selected with biomartr::getMarts()
                                dataset    = "athaliana_eg_gene",               # datasets were selected with biomartr::getDatasets()
                                attributes = attributes_to_retrieve,            # attributes were selected with biomartr::getAttributes()
                                filters = 	"ensembl_gene_id" )# query key


##########################
# 03. ORA Cluster Profiler
##########################


##############################
# 03. ORA with interproscan GO
##############################

interpro <- read.delim("00.tutorial/interproscan_results.tsv", header = F, stringsAsFactors = F, row.names = NULL)

new_colnames <- c("protein_id",
                  "md5",
                  "seq_len",
                  "analysis",
                  "signature_accession",
                  "signature_description",
                  "start",
                  "stop",
                  "score",
                  "status",
                  "date",
                  "interpro_accession",
                  "interpro_description",
                  "go")
colnames(interpro) <- new_colnames

interpro_go <- 
  interpro %>% 
  dplyr::select(protein_id, go) %>% 
  dplyr::filter(go != "-") %>% 
  dplyr::filter(go != "")
tail(interpro_go, n = 10)

splitted_interpro_go <- cSplit(indt = interpro_go, splitCols = "go", sep = "|", direction = "long")
dedup_go  <- splitted_interpro_go %>% distinct()
tail(dedup_go)

# change column name
# remove the .1, .2 at the end of the gene
dedup_go <- dedup_go %>% 
  dplyr::rename("genes" = "protein_id") %>% 
  mutate(genes = substr(genes, 1, 9))
tail(dedup_go)

write.csv(x = dedup_go, file = "00.tutorial/gene_ontologies_all_genes.csv", row.names = F)

# subset diff genes
diff_genes <- read_delim(file = "00.tutorial/differential_genes.tsv", delim = "\t") 
  

diff_genes_go <- inner_join(x = dedup_go, y = diff_genes)
head(diff_genes_go)

diff_genes_go %>% 
  select(genes, go) %>% 
  write.csv(file = "00.tutorial/gene_ontologies_diff_genes.csv", row.names = F)

#################################################
# KEGG analysis using custom KEGG classification
#################################################

transcriptome <- read.table("00.tutorial/Araport11_genes_ko.txt", fill = T, sep = "\t", row.names = NULL, header = F)
differential_genes <- read.table("00.tutorial/differential_genes.tsv", sep = "\t", header = T)

transcriptome[,1] <- gsub("\\..*","",transcriptome[,1])
transcriptome <- distinct(transcriptome)

KITotal <- as.data.frame(table(transcriptome$V2))
KITotal <- KITotal[-1,]
bgTotal <- sum(KITotal[,2])

transcriptome_diff_filt <- transcriptome[which(transcriptome$V1 %in% differential_genes$genes),]
KIquery <- as.data.frame(table(transcriptome_diff_filt$V2))
KIquery <- KIquery[-1,]
queryTotal <- sum(KIquery[,2])


tableForPhyper <- merge(KIquery, KITotal, by = "Var1")
colnames(tableForPhyper)[1:3] <- c("KOID","queryItem","bgItem")
tableForPhyper$queryTotal <- queryTotal
tableForPhyper$bgTotal <- bgTotal

res <- tableForPhyper %>%                                  
  rowwise() %>% 
  mutate(p_val = phyper(q = queryItem, m = bgItem, n = bgTotal-queryTotal, k = queryTotal))
