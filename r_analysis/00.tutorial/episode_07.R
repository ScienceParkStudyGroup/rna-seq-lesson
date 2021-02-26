suppressPackageStartupMessages(library("biomartr"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("org.At.tair.db"))
suppressPackageStartupMessages(library("biomaRt"))  # only use to remove cache bug


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
#head(result_BM)  

##########################
# 03. ORA Cluster Profiler
##########################


