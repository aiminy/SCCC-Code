Re.Jun<-read.csv("GeneWise_jscs3_all_with_anno_2_24_2016.csv")
Re.DE<-read.csv("DE_cheng_output_sample1_DE.csv")

head(Re.Jun)
head(Re.DE)

Re.Jun[,3]

mart = useMart('ensembl')
listDatasets(mart)

ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
listFilters(ensembl)


#geneMapsL<-array()

#for(i in 1:length(Re.DE[,1])){
geneMaps<-getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters="mgi_symbol",values=Re.DE[,1],mart=ensembl)

geneMaps<-getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters="mgi_symbol",values=Re.DE[,1],mart=ensembl)

dim(geneMaps)
length(Re.DE[,1])
head(geneMaps)

head(geneMaps)
dim(geneMaps)

colnames(Re.DE)[1]="mgi_symbol"

Re.DE2<-merge(geneMaps,Re.DE,by="mgi_symbol")

head(Re.DE2)
colnames(Re.DE2)[2]="geneID"


Re.DE.Jun<-merge(Re.DE2,Re.Jun,by="geneID")

dim(Re.DE.Jun)



head(Re.DE.Jun)
colnames(Re.DE.Jun)

plot(Re.DE.Jun[,8],Re.DE.Jun[,23])

hist(Re.DE.Jun[,8])
hist(Re.DE.Jun[,23])


which(is.na(geneMaps[,1]))
length(unique(Re.DE[,1]))

cbind(geneMaps,Re.DE[,1])

#if(length(geneMaps)==0){
#  geneMaps="NA"
#}
  
#geneMapsL<-c(geneMapsL,geneMaps)

#}

geneMaps.2<-getBM(attributes=c("mgi_symbol"), filters="ensembl_gene_id",values=geneMaps,mart=ensembl)

geneMaps.3<-cbind2(geneMaps,geneMaps.2)


Re.DE.with.anno<-cbind(geneMaps,Re.DE)


mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
getBM(attributes=c("ensembl_gene_id"), filters="ensembl_gene_id",values=allGenes.hallmark,mart=mart)

getBM(attributes=c("ensembl_gene_id"),values=allGenes.hallmark,mart=mart)



function(data.gene){
  gene.anno<-as.character(paste(getBM(attributes=c("mgi_symbol","description"), filters="ensembl_gene_id", 
                                      values=data.gene, 
                                      mart=mart)$mgi_symbol,
                                getBM(attributes=c("mgi_symbol","description"), filters="ensembl_gene_id", 
                                      values=data.gene, 
                                      mart=mart)$description,sep="->"))
  return(gene.anno)
}


save.image(file="/media/H_driver/2015/Nimer_Cheng/Data_set_two_methods.RData")
savehistory(file="/media/H_driver/2015/Nimer_Cheng/Data_set_two_methods.Rhistory")

