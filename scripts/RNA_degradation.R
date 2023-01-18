library(dplyr)
path<-"results/study9/htseq_count_exon"
fileNames<-dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   
exon.df.list <- lapply(filePath, function(x){
  read.delim(x,header = F)})

CalDegradategenes <- function(exon.df) {
  # for each sample, calculate the # of degraded genes / total genes
  uni_id = unique(exon.df$geneid)
  
  
  degrade.id <- 0
  for (i in 1:length(uni_id)) {
    target.id <-  uni_id[i]
    
    df1 <- exon.df %>% filter(geneid == target.id)
    if (sum(df1[df1$geneid ==  target.id, ]["count"] > 0) == 1) {
      df2 <- df1 %>% filter(exon_number == max(exon_number)) %>%
        filter(count != 0)
      if (nrow(df2) > 0) {
        degrade.id = degrade.id + 1
        #print(target.id)
      }
    }
  }
  return(degrade.id)
}


RNA_degradation_genes <- data.frame()

for (k in 1:length(exon.df.list)) {
  exon.df = exon.df.list[[k]]
  SRR.name <- gsub("\\.count","",fileNames[k])
  
  colnames(exon.df) <-
    c("geneid:exon_number", "gene_name", "exon_number", "count")
  exon.df$`geneid:exon_number` <-
    gsub("\\:\\d*", "", exon.df$`geneid:exon_number`)
  colnames(exon.df)[1] <- "geneid"
  print(fileNames[k])
  
  tmp.number <- CalDegradategenes(exon.df)
  
  RNA_degradation_genes[k, 1] <- SRR.name
  RNA_degradation_genes[k, 2] <- tmp.number
  colnames(RNA_degradation_genes) <- c("Run ID", "number_RNA_degration_genes")
}

write.csv(RNA_degradation_genes,file="results/study9/RNA_degradation/degrate_genes.csv")
