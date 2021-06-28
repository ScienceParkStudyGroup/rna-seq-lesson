library("org.At.tair.db") # version 3.10.0

columns(org.At.tair.db)
#[1] "ARACYC"       "ARACYCENZYME" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
#[7] "GENENAME"     "GO"           "GOALL"        "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
#[13] "PMID"         "REFSEQ"       "SYMBOL"       "TAIR"        

length(keys(org.At.tair.db, keytype = 'TAIR')) # 27416

length(keys(org.At.tair.db, keytype = 'GO')) # 4837