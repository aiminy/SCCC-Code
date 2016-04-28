#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

biocLite("mogene20sttranscriptcluster.db")
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

biocLite("pd.hg.u133.plus.2")
library(pd.hg.u133.plus.2)
biocLite("ggplot2")
library(ggplot2)
biocLite("reshape")
biocLite("plot3D")
biocLite("gplots")
biocLite("ggdendro")
biocLite("RColorBrewer")
biocLite("frma")

library("reshape")
library("plot3D")
library("gplots")
library(ggdendro)
library(RColorBrewer)
library(RColorBrewer)
library(frma)
biocLite("hgu133plus2frmavecs")
#R CMD INSTALL -l R/x86_64-pc-linux-gnu-library/3.2/ /media/H_driver/Aimin_project/hgu133plus2frmavecs_1.5.0.tar.gz
library(hgu133plus2frmavecs)
biocLite("Homo.sapiens")
library("Homo.sapiens")

# Read in the CEL files in the directory, then normalize the data
read.celfiles(list.celfiles("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM"))
affybatch <- read.celfiles(list.celfiles("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM"))
eset <- rma(affybatch)

Function2ReadCEL_or_cel<-function(input_dir){
  affybatch <- read.celfiles(list.celfiles(input_dir))
  eset <- rma(affybatch)
  re<-list(aff_data=affybatch,ed_norma=eset)
  return(re)
}

Re.set1<-Function2ReadCEL_or_cel("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")


setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
CEL.filenames.data.set1=dir(pattern="CEL","/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy(filenames = CEL.filenames.data.set1)
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)
ed.normalized.set1<- rma(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)
ed.normalized.set2<- rma(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)
ed.normalized.set3<- rma(data.set3)

data.set1.normalized<-ed.normalized.set1
data.set2.normalized<-ed.normalized.set2
data.set3.normalized<-ed.normalized.set3

ll(dim=T)

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")


cel.file.all<-c(unique(colnames(data.set1.normalized)),unique(colnames(data.set2.normalized)),unique(colnames(data.set3.normalized)))

length(cel.file.all)
cel.file.all

cel.file.all.2<-rbind(cbind(unique(colnames(data.set1.normalized)),rep(1,length(unique(colnames(data.set1.normalized))))),
cbind(unique(colnames(data.set2.normalized)),rep(2,length(unique(colnames(data.set2.normalized))))),
cbind(unique(colnames(data.set3.normalized)),rep(3,length(unique(colnames(data.set3.normalized))))))

cel.file.cancer.data.set2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,8])
cel.file.cancer.data.set2.2<-as.character(data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)][,1])
cel.file.cancer.data.set2.2.2<-cel.file.cancer.data.set2.2[-which(cel.file.cancer.data.set2.2=="")]

data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.cel.mapping.3<-data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),]

data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",as.character(data.set2.cancer.GSM.cel.mapping.3[,2]))
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))

cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])

cel.file.cancer.set1.set2<-rbind(cbind(cel.file.cancer.data.set1,rep(1,length(cel.file.cancer.data.set1))),
cbind(cel.file.cancer.data.set2,rep(2,length(cel.file.cancer.data.set2))),
cbind(cel.file.cancer.data.set2.2.2,rep(2,length(cel.file.cancer.data.set2.2.2))))

cel.file.cancer.set1.set2.set3<-rbind(cel.file.cancer.set1.set2,cel.file.all.2[which(cel.file.all.2[,2]==3),])

cel.file.cancer.set1.set2.set3.name<-gsub(".CEL","",cel.file.cancer.set1.set2.set3[,1])
cel.file.cancer.set1.set2.set3.name.1<-gsub(".cel","",cel.file.cancer.set1.set2.set3.name)
cel.file.cancer.set1.set2.set3.name.2<-gsub("-","_",cel.file.cancer.set1.set2.set3.name.1)

data.set123.normalized<-cbind(data.set1.normalized,data.set2.normalized,data.set3.normalized)

dim(data.set123.normalized)
dim(data.set123.normalized[,which(colnames(data.set123.normalized) %in% cel.file.cancer.set1.set2.set3[,1])])
colnames(data.set123.normalized)[grep(54,colnames(data.set123.normalized))]

original.cel.file.names<-colnames(data.set123.normalized)
original.cel.file.names.1<-gsub("X","",original.cel.file.names)
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.4<-gsub("\\.","_",original.cel.file.names.3)

index.4.cancer.sample<-which(original.cel.file.names.4 %in% cel.file.cancer.set1.set2.set3.name.2)
original.cel.file.name.with.mapped.names<-cbind(original.cel.file.names[index.4.cancer.sample],original.cel.file.names.4[index.4.cancer.sample])
#index.4.no.cancer.sample.1<-cbind(original.cel.file.names[-index.4.cancer.sample],original.cel.file.names.4[-index.4.cancer.sample])
setdiff(cel.file.cancer.set1.set2.set3.name.2,original.cel.file.names.4[index.4.cancer.sample])
which(original.cel.file.names.4[-index.4.cancer.sample] %in% cel.file.cancer.set1.set2.set3.name.2)

cancer.data.set123.normalized<-data.set123.normalized[,index.4.cancer.sample]

cancer.cel.file.name.72<-colnames(cancer.data.set123.normalized)
cancer.cel.file.name.72.1<-gsub("X","",cancer.cel.file.name.72)
cancer.cel.file.name.72.2<-gsub(".CEL","",cancer.cel.file.name.72.1)
cancer.cel.file.name.72.3<-gsub(".cel","",cancer.cel.file.name.72.2)
cancer.cel.file.name.72.4<-gsub("\\.","_",cancer.cel.file.name.72.3)
cancer.cel.file.name.72.5<-cancer.cel.file.name.72.4

data.set2.cancer.GSM.cel.mapping.44<-as.character(data.set2.cancer.GSM.cel.mapping.4[which(data.set2.cancer.GSM.cel.mapping.4[,1] %in% cancer.cel.file.name.72.4),3])
cancer.cel.file.name.72.5[which(cancer.cel.file.name.72.5 %in% data.set2.cancer.GSM.cel.mapping.4[,1])]<-data.set2.cancer.GSM.cel.mapping.44
cancer.cel.file.name.72.6<-cbind(cancer.cel.file.name.72,cancer.cel.file.name.72.5)
cancer.cel.file.name.72.7<-c(sapply(strsplit(cancer.cel.file.name.72.6[1:29,2],"_"),"[[",1),
sapply(strsplit(cancer.cel.file.name.72.6[30:72,2],"_"),"[[",4))
cancer.cel.file.name.72.8<-cbind(cancer.cel.file.name.72.6,cancer.cel.file.name.72.7)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.normalized.st1<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.1)]
cancer.data.set123.normalized.st2<-as.data.frame(cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.2)])
colnames(cancer.data.set123.normalized.st2)<-colnames(cancer.data.set123.normalized)[which(cancer.cel.file.name.72.8[,3] %in% subtype.2)]
cancer.data.set123.normalized.st3<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.3)]
cancer.data.set123.normalized.st4<-cancer.data.set123.normalized[,which(cancer.cel.file.name.72.8[,3] %in% subtype.4)]

#n.st<-length(colnames(cancer.data.set123.normalized.st1))

# #Use subtype1,2,3,4
# cel.file.sample.infor<-as.data.frame(rbind(
# cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
# cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
# cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
# cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
# ))
# colnames(cel.file.sample.infor)=c("filename","subtype")

#Use subtype 1,3,4
cel.file.sample.infor.no.2<-as.data.frame(rbind(
  cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
  cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
  cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))

#Use frma method to do normalization  
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

cel.file.sample.infor.no.3<-gsub("X","",cel.file.sample.infor.no.2[,1])
cel.file.sample.infor.no.4<-gsub(".CEL","",cel.file.sample.infor.no.3)
cel.file.sample.infor.no.5<-gsub(".cel","",cel.file.sample.infor.no.4)
cel.file.sample.infor.no.6<-gsub("\\.","_",cel.file.sample.infor.no.5)
cel.file.sample.infor.no.7<-cbind(cel.file.sample.infor.no.2,cel.file.sample.infor.no.6)
colnames(cel.file.sample.infor.no.7)=c("filename","subtype","match_key")

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.sample.infor.no.7[,3])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

setdiff(cel.file.sample.infor.no.3[,3],samp.frma.2[,2])
samp.frma.1[grep(1443,samp.frma.1[,2]),]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]

dim(ed.set123.frma.cancer)
samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("match_key","filename_frma")
cel.file.sample.infor.no.8<-merge(cel.file.sample.infor.no.7,samp.frma.8,by="match_key")

head(samp.frma.7)
dim(ed.set123.frma.cancer)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")
ndata<-ed.set123.frma.cancer
geneSymbol=GeneSym.all
tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

rownames(tempdata.byGSym.2) = NULL
data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,max)
)

data.byGSym.index = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,which.max)
)

#class(data.byGSym)

#which(is.na(data.byGSym$Symbol))
rownames(data.byGSym)=data.byGSym$Symbol
data.byGSym.2<-data.byGSym[,-61]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp
SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma.pdf","PCA_allsample_frma.pdf","/media/H_driver/2015/Sophia/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

fit.st134.gene.frma <- lmFit(data.byGSym.2, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=60000)
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=60000)


OutPut2HtmlTable<-function(TopTableSt41.gene.frma,out_dir,out_file_name,out_title){

hgnc.gene.symbol.ENTREZID<-select(Homo.sapiens, keys=rownames(TopTableSt41.gene.frma),columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")

TopTableSt41.gene.frma.2<-data.frame(rownames(TopTableSt41.gene.frma),TopTableSt41.gene.frma)
colnames(TopTableSt41.gene.frma.2)[1]="SYMBOL"

TopTableSt41.gene.frma.4<-merge(TopTableSt41.gene.frma.2,hgnc.gene.symbol.ENTREZID,by="SYMBOL",sort=FALSE)
htmlpage(list(TopTableSt41.gene.frma.4$ENTREZID),filename=paste0(out_dir,out_file_name), 
         title=out_title, 
         othernames=TopTableSt41.gene.frma.4[,-8], 
         table.head=c("ENTREZID",colnames(TopTableSt41.gene.frma.3[,-8])), 
         table.center=TRUE, digits=6)
}

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/","St41-5-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/","St34-5-frma.html","st3-vs-st4_HTML_report")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
