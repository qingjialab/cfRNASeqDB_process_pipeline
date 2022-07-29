#configfile: "envs/config.yaml"
#STUDY = config["study5"]["name"]
#RESULT_DIR = "results/" + STUDY + "/"

hg_symbols<-read.csv("rRNA/hg_symbols.csv")
path<-"results/study9/htseq_count"
fileNames<-dir(path)
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})   
data <- lapply(filePath, function(x){
  read.table(x)}) 
rRNA_rate<-data.frame()
for (i in 1:50 ) {
  data[[i]]<-data.frame(data[[i]])[-61542:-61546,]
  colnames(data[[i]])<-c("gene_id","count")
  a<-merge(data[[i]],hg_symbols,by="gene_id")
  gene_RP_MT<-a[grep("ribosomal",a$description),]
  rate<-colSums(data.frame(gene_RP_MT$count))/colSums(data.frame(data[[i]]$count))
  name<-names(data[i])
  name<-substring(name,1,nchar(name)-6)
  rRNA_rate[i,1]<-name
  rRNA_rate[i,2]<-rate
  colnames(rRNA_rate)<-c("sample","rRNA rate")
}
write.csv(rRNA_rate,file="results/study9/rRNA_percentage/rRNA_rate.csv")
