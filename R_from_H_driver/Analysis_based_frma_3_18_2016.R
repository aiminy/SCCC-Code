#load the affy library
.libPaths()
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(plyr)

#biocLite("mogene20sttranscriptcluster.db")
#library(mogene20stprobeset.db)
#library(mogene20sttranscriptcluster.db)

#biocLite("hugene10sttranscriptcluster.db")
#library(hugene10sttranscriptcluster.db)

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
library(hgu133plus2.db)

# Read the CEL files in set1,2 and set3(cancer)
setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
ed.raw.set1 <- exprs(data.set1)
samp.set1 <- sampleNames(data.set1)
probes.set1 <- featureNames(data.set1)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set2_corrected")
data.set2 <- ReadAffy()
ed.raw.set2 <- exprs(data.set2)
samp.set2 <- sampleNames(data.set2)
probes.set2 <- featureNames(data.set2)

setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set3_HGSCancer_only")
data.set3 <- ReadAffy()
ed.raw.set3 <- exprs(data.set3)
samp.set3 <- sampleNames(data.set3)
probes.set3 <- featureNames(data.set3)

ll(dim=T)

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Set2_header_missing/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")

data.ER.PR.sample.info.new<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx".sheet=2)

data.ER.PR.sample.info.sheet1<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx",sheet=1)
data.ER.PR.sample.info.sheet2<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx",sheet=2)

original.cel.file.all<-rbind(cbind(unique(samp.set1),rep(1,length(unique(samp.set1)))),
cbind(unique(samp.set2),rep(2,length(unique(samp.set2)))),
cbind(unique(samp.set3),rep(3,length(unique(samp.set3)))))

original.cel.file.names.1<-gsub(" - ","_",original.cel.file.all[,1])
original.cel.file.names.2<-gsub(".CEL","",original.cel.file.names.1)
original.cel.file.names.3<-gsub(".cel","",original.cel.file.names.2)
original.cel.file.names.key<-cbind(original.cel.file.all,gsub("-","_",original.cel.file.names.3))
colnames(original.cel.file.names.key)=c("cel_file_name","set","key")

#Extract all cel files belonging to the cancer samples in dat set2
data.set2.cancer.GSM.cel.mapping<-data.set2.sample.info[which(data.set2.sample.info[,13]=="cancer"),c(1:12,13)]
data.set2.cancer.GSM.cel.mapping.2<-data.set2.cancer.GSM.cel.mapping[,c(1,8)]
data.set2.cancer.GSM.part<-as.character(data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),1])
data.set2.cancer.GSM.part.cel.name<-as.character(data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])!=""),2])

ReformatCelName<-function(data.set2.cancer.GSM.cel.mapping.3){
data.set2.cancer.GSM.cel.mapping.3.2<-gsub(".cel","",data.set2.cancer.GSM.cel.mapping.3)
data.set2.cancer.GSM.cel.mapping.3.2.2<-gsub(".CEL","",as.character(data.set2.cancer.GSM.cel.mapping.3.2))
data.set2.cancer.GSM.cel.mapping.3.2.2.2<-gsub("-","_",as.character(data.set2.cancer.GSM.cel.mapping.3.2.2))
data.set2.cancer.GSM.cel.mapping.4<-cbind(data.set2.cancer.GSM.cel.mapping.3,data.set2.cancer.GSM.cel.mapping.3.2.2.2,sapply(strsplit(data.set2.cancer.GSM.cel.mapping.3.2.2.2,"_"),"[[",1))
return(data.set2.cancer.GSM.cel.mapping.4)
}

data.set2.cancer.GSM.cel.mapping.combine<-cbind(data.set2.cancer.GSM.part,ReformatCelName(data.set2.cancer.GSM.part.cel.name))
GSM.based.cel.index<-which(original.cel.file.names.key[,3] %in% data.set2.cancer.GSM.cel.mapping.combine[,1])
cel.key.based.index<-which(original.cel.file.names.key[,3] %in% data.set2.cancer.GSM.cel.mapping.combine[,3])
com.GSM.cel.index<-intersect(GSM.based.cel.index,cel.key.based.index)
diff.GSM.cel.index<-setdiff(GSM.based.cel.index,cel.key.based.index)
diff.cel.GSM.index<-setdiff(cel.key.based.index,GSM.based.cel.index)
original.cel.file.names.key.rm.dup<-original.cel.file.names.key[-diff.cel.GSM.index,]

original.cel.file.names.set2.cancer.index<-which(original.cel.file.names.key.rm.dup[,3] %in% data.set2.cancer.GSM.cel.mapping.combine[,1])
original.cel.file.names.set2.cancer<-original.cel.file.names.key.rm.dup[original.cel.file.names.set2.cancer.index,]
data.set2.cancer.cel.part.name<-as.character(data.set2.cancer.GSM.cel.mapping.2[which(as.character(data.set2.cancer.GSM.cel.mapping.2[,1])==""),2])
data.set2.cancer.cel.part.name.combine<-cbind(data.set2.cancer.cel.part.name,ReformatCelName(data.set2.cancer.cel.part.name))
original.cel.file.names.set2.cancer2.index<-which(original.cel.file.names.key.rm.dup[,3] %in% data.set2.cancer.cel.part.name.combine[,3])
original.cel.file.names.set2.cancer2<-original.cel.file.names.key.rm.dup[original.cel.file.names.set2.cancer2.index,]
original.cel.file.names.set2.cancer.all<-rbind(original.cel.file.names.set2.cancer,original.cel.file.names.set2.cancer2)
colnames(original.cel.file.names.set2.cancer.all)=c("cel_file_name","set","key")
data.set2.cancer.GSM.cel.key<-rbind(data.set2.cancer.GSM.cel.mapping.combine[,c(1,4)],data.set2.cancer.cel.part.name.combine[,c(3,4)])
colnames(data.set2.cancer.GSM.cel.key)=c("key","symbol")
original.cel.file.names.set2.cancer.all.2<-merge(original.cel.file.names.set2.cancer.all,data.set2.cancer.GSM.cel.key,by="key",sort=FALSE)

#Extract all cel files belonging to the cancer samples in dat set1
cel.file.cancer.data.set1<-as.character(data.set1.sample.info[which(data.set1.sample.info[,4]=="cancer"),2])
data.set1.cancer.celname.combine<-cbind(cel.file.cancer.data.set1,ReformatCelName(cel.file.cancer.data.set1))
data.set1.cancer.celname.combine.2<-data.set1.cancer.celname.combine[,3:4]
colnames(data.set1.cancer.celname.combine.2)=c("key","symbol")
original.cel.file.names.set1.cancer<-merge(original.cel.file.names.key,data.set1.cancer.celname.combine.2,by="key",sort=FALSE)

#Generate the symbol for all cel files belonging to the cancer samples in dat set3
original.cel.file.names.key.data.set3<-original.cel.file.names.key[which(original.cel.file.names.key[,2]==3),]
original.cel.file.names.key.data.set3.2<-cbind(original.cel.file.names.key.data.set3,sapply(strsplit(original.cel.file.names.key.data.set3[,3],"_"),"[[",4))
colnames(original.cel.file.names.key.data.set3.2)[4]="symbol"
original.cel.file.names.key.data.set3.3<-original.cel.file.names.key.data.set3.2[,c(3,1,2,4)]

#Combine the cancer samples in date set1,2,3
cel.file.cancer.set1.set2.set3<-rbind(original.cel.file.names.set1.cancer,original.cel.file.names.set2.cancer.all.2,original.cel.file.names.key.data.set3.3)

#Classify cancer patients to 4 subtypes

subtype.1<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),1])
subtype.2<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),1])
subtype.3<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),1])
subtype.4<-as.character(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),1])

cancer.data.set123.st1<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.1),]
cancer.data.set123.st2<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.2),]
cancer.data.set123.st3<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.3),]
cancer.data.set123.st4<-cel.file.cancer.set1.set2.set3[which(cel.file.cancer.set1.set2.set3[,4] %in% subtype.4),]

cancer.data.set123.st1.subtype<-cbind(cancer.data.set123.st1,rep("st1",dim(cancer.data.set123.st1)[1]))
colnames(cancer.data.set123.st1.subtype)[5]="subtype"
cancer.data.set123.st2.subtype<-cbind(cancer.data.set123.st2,rep("st2",dim(cancer.data.set123.st2)[1]))
colnames(cancer.data.set123.st2.subtype)[5]="subtype"
cancer.data.set123.st3.subtype<-cbind(cancer.data.set123.st3,rep("st3",dim(cancer.data.set123.st3)[1]))
colnames(cancer.data.set123.st3.subtype)[5]="subtype"
cancer.data.set123.st4.subtype<-cbind(cancer.data.set123.st4,rep("st4",dim(cancer.data.set123.st4)[1]))
colnames(cancer.data.set123.st4.subtype)[5]="subtype"

#Use subtype 1,2,3,4
cel.file.name.key.set.symbol.subtype<-rbind(cancer.data.set123.st1.subtype,cancer.data.set123.st2.subtype,cancer.data.set123.st3.subtype,cancer.data.set123.st4.subtype)

#Use subtype 1,3,4
cel.file.name.key.set.symbol.st134<-cel.file.name.key.set.symbol.subtype[which(cel.file.name.key.set.symbol.subtype[,5]!="st2"),]

#Use frma method to do normalization
setwd("/media/H_driver/2015/Sophia/Cel_file_frma/")
data.set123.raw <- ReadAffy()
data.set123.frma<-frma(data.set123.raw)
ed.set123.frma <- exprs(data.set123.frma)
samp.frma <- sampleNames(data.set123.frma)
probes.frma <- featureNames(data.set123.frma)

samp.frma.1<-gsub(" - ","_",samp.frma)
samp.frma.2<-gsub(".CEL","",samp.frma.1)
samp.frma.3<-gsub(".cel","",samp.frma.2)
samp.frma.4<-gsub("-","_",samp.frma.3)
samp.frma.5<-cbind(samp.frma,samp.frma.4)
samp.frma.cancer.index<-which(samp.frma.5[,2] %in% cel.file.name.key.set.symbol.st134[,1])
samp.frma.6<-samp.frma.5[samp.frma.cancer.index,]

data.set123.frma.cancer<-data.set123.frma[,samp.frma.cancer.index]
ed.set123.frma.cancer<-ed.set123.frma[,samp.frma.cancer.index]
dim(ed.set123.frma.cancer)

samp.frma.7<-cbind(samp.frma.6,colnames(ed.set123.frma.cancer))
samp.frma.8<-samp.frma.7[,2:3]
colnames(samp.frma.8)<-c("key","filename_frma")

cel.file.sample.infor.no.8<-merge(cel.file.name.key.set.symbol.st134,samp.frma.8,by="key",sort=FALSE)

GeneSym.all <- getSYMBOL(rownames(ed.set123.frma.cancer), "hgu133plus2.db")


#<-cbind(names(GeneSym.all),GeneSym.all)

# #Collapse probes into genes
ndata<-ed.set123.frma.cancer

geneSymbol=GeneSym.all
tempdata.byGSym =data.frame(ndata, Symbol = geneSymbol)

colnames(tempdata.byGSym)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

# test.sample<-ed.set123.frma.cancer[which(GeneSym.all=="ERBB4"),which(colnames(ed.set123.frma.cancer) %in% cel.file.sample.infor.no.8[which(cel.file.sample.infor.no.8[,5]=="st4"),6])]
# 
# test.sample.2<-data.frame(test.sample,Symbol="ERBB4")
# test.sample.3 = test.sample.2
# test.sample.3[,1:5] = apply(test.sample.3[,1:5],2,as.numeric)

data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h[,1:(dim(tempdata.byGSym.2)[2]-1)],2,max)
)


rownames(data.byGSym)=data.byGSym$Symbol
colnames(data.byGSym)

data.byGSym.2<-data.byGSym[,-1]
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
colnames(temp) = colnames(ndata)

data.byGSym.2 = temp
colnames(data.byGSym.2)
#data.byGSym.2 = ed.set123.frma.cancer

SampleType<-cel.file.sample.infor.no.8$subtype
heatmap_wPCA(data.byGSym.2, "heatmap_allsample_frma_all_probes.pdf","PCA_allsample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/", SampleType)

f.st134.frma <- factor(cel.file.sample.infor.no.8$subtype)
design.st134.frma <- model.matrix(~0+f.st134.frma)
colnames(design.st134.frma) <- levels(f.st134.frma)

reorder.index<-match(cel.file.sample.infor.no.8[,6],colnames(data.byGSym.2))

data.byGSym.3<-data.byGSym.2[,reorder.index]

#colnames(data.byGSym.3)

fit.st134.gene.frma <- lmFit(data.byGSym.3, design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
fit2.st134.frma  <- contrasts.fit(fit.st134.gene.frma, cont.matrix.st134.frma)
fit2.st134.frma  <- eBayes(fit2.st134.frma)

TopTableSt34.gene.frma<-topTable(fit2.st134.frma,coef=1,n=dim(data.byGSym.2)[1])
TopTableSt41.gene.frma<-topTable(fit2.st134.frma,coef=2,n=dim(data.byGSym.2)[1])

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

OutPut2HtmlTable(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-8-frma.html","st4-vs-st1_HTML_report")
OutPut2HtmlTable(TopTableSt34.gene.frma,"/media/H_driver/2015/Sophia/Results/","St34-8-frma.html","st3-vs-st4_HTML_report")

#Generate the rank file for GSEA
GenerateRankFile4GSEA<-function(TopTableSt41.gene.frma,out_dir,out_file_name,up_down_top){

  FC.sign=sign(TopTableSt41.gene.frma[,1])
  p.inve=1/TopTableSt41.gene.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve

  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.gene.frma,Gene.Rank.Score)

  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]

  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("GeneName","Score")

  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))

  write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)

  return(TopTableSt41.gene.frma.GSEA.sorted.by.score)
}

GenerateRankFile4GSEA(TopTableSt41.gene.frma,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA_up_down_300.rnk",300)

#write.csv(TopTableSt41.gene.frma, file="/media/H_driver/2015/Sophia/Results/St41.csv",quote =FALSE)


##Use all probes for DE 
reorder.index.all.probes<-match(cel.file.sample.infor.no.8[,6],colnames(ed.set123.frma.cancer))
ed.set123.frma.cancer.reorder<-ed.set123.frma.cancer[,reorder.index.all.probes]

check.match<-cbind(colnames(ed.set123.frma.cancer.reorder),design.st134.frma,cel.file.sample.infor.no.8)

fit.st134.all.probes.frma <- lmFit(ed.set123.frma.cancer.reorder, design.st134.frma)
#cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134.frma)
cont.matrix.st134.frma <- makeContrasts(st34="st3-st4",st41="st4-st1",st13="st1-st3",levels=design.st134.frma)

fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
fit2.st134.all.probes.frma  <- eBayes(fit2.st134.all.probes.frma)

TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(ed.set123.frma.cancer.reorder)[1])
TopTableSt41.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=2,n=dim(ed.set123.frma.cancer.reorder)[1])
TopTableSt13.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=3,n=dim(ed.set123.frma.cancer.reorder)[1])

OutPut2HtmlTableAllProbes<-function(TopTableSt41.all.probes.frma,out_dir,out_file_name,out_title){
  
  
  GeneSym.all <- as.data.frame(getSYMBOL(rownames(TopTableSt41.all.probes.frma), "hgu133plus2.db"))
  
  #probes.GeneSym.all<-cbind(rownames(TopTableSt41.all.probes.frma),GeneSym.all)
  
  #rownames(probes.GeneSym.all)
  #hgnc.gene.symbol.ENTREZID<-select(Homo.sapiens, keys=rownames(TopTableSt41.gene.frma),columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
  
  
  TopTableSt41.all.probes.frma.2<-merge(GeneSym.all,TopTableSt41.all.probes.frma,by=0,sort=FALSE)
  rownames(TopTableSt41.all.probes.frma.2)=TopTableSt41.all.probes.frma.2[,1]
  
  colnames(TopTableSt41.all.probes.frma.2)[1]="Probe_ID" 
  colnames(TopTableSt41.all.probes.frma.2)[2]="SYMBOL"
  
 TopTableSt41.all.probes.frma.3<-data.frame(TopTableSt41.all.probes.frma.2,Index=seq(1,dim(TopTableSt41.all.probes.frma.2)[1]))
  
  #colnames(TopTableSt41.gene.frma.2)[1]="SYMBOL"
  ll <- getEG(rownames(TopTableSt41.all.probes.frma.3),"hgu133plus2.db")
  
  #TopTableSt41.gene.frma.4<-merge(TopTableSt41.gene.frma.2,hgnc.gene.symbol.ENTREZID,by="SYMBOL",sort=FALSE)

    
  htmlpage(list(ll),filename=paste0(out_dir,out_file_name),
           title=out_title,
           othernames=TopTableSt41.all.probes.frma.3,
           table.head=c("Locus_ID",colnames(TopTableSt41.all.probes.frma.3)),
           table.center=TRUE, digits=6)
  
  return(TopTableSt41.all.probes.frma.3)
}

RE41<-OutPut2HtmlTableAllProbes(TopTableSt41.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St41-all-probes-frma-3.html","st4-vs-st1_HTML_report")
RE34<-OutPut2HtmlTableAllProbes(TopTableSt34.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St34-all-probes-frma-3.html","st3-vs-st4_HTML_report")
RE13<-OutPut2HtmlTableAllProbes(TopTableSt13.all.probes.frma,"/media/H_driver/2015/Sophia/Results/","St13-all-probes-frma.html","st1-vs-st3_HTML_report")

#head(RE41)

GenerateRank4AllProbes<-function(TopTableSt41.all.probes.frma,up_down_top){
  
  FC.sign=sign(TopTableSt41.all.probes.frma[,1])
  p.inve=1/TopTableSt41.all.probes.frma[,4]
  Gene.Rank.Score=FC.sign*p.inve
  
  TopTableSt41.gene.frma.GSEA<-cbind(TopTableSt41.all.probes.frma,Gene.Rank.Score)
  
  TopTableSt41.gene.frma.GSEA.sorted.by.score<-TopTableSt41.gene.frma.GSEA[order(TopTableSt41.gene.frma.GSEA[,7],decreasing = TRUE),]
  
  Re.GSEA<-cbind(rownames(TopTableSt41.gene.frma.GSEA.sorted.by.score),TopTableSt41.gene.frma.GSEA.sorted.by.score[,7])
  colnames(Re.GSEA)=c("Probe","Score")
  
  Re.GSEA.2<-rbind(head(Re.GSEA,up_down_top),tail(Re.GSEA,up_down_top))
  
  #write.table(Re.GSEA.2, file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = FALSE)
  
  return(Re.GSEA.2)
}


TopTableSt41.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt41.all.probes.frma,200)
TopTableSt34.most.DE.probes.frma<-GenerateRank4AllProbes(TopTableSt34.all.probes.frma,200)

ed.set123.frma.cancer.reorder.st1<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st1==1),1])]
ed.set123.frma.cancer.reorder.st2<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st2==1),1])]
ed.set123.frma.cancer.reorder.st3<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st3==1),1])]
ed.set123.frma.cancer.reorder.st4<-ed.set123.frma.cancer.reorder[,which(colnames(ed.set123.frma.cancer.reorder) %in% check.match[which(check.match$st4==1),1])]

dim(ed.set123.frma.cancer.reorder.st1)
dim(ed.set123.frma.cancer.reorder.st2)
dim(ed.set123.frma.cancer.reorder.st3)
dim(ed.set123.frma.cancer.reorder.st4)

ed.set123.frma.cancer.reorder.st41<-cbind(ed.set123.frma.cancer.reorder.st4,ed.set123.frma.cancer.reorder.st1)
ed.set123.frma.cancer.reorder.st34<-cbind(ed.set123.frma.cancer.reorder.st3,ed.set123.frma.cancer.reorder.st4)

dim(ed.set123.frma.cancer.reorder.st41)
dim(ed.set123.frma.cancer.reorder.st34)

ed.set123.frma.cancer.reorder.st41.MDE.100<-ed.set123.frma.cancer.reorder.st41[which(rownames(ed.set123.frma.cancer.reorder.st41) %in%  TopTableSt41.most.DE.probes.frma[,1]),]
ed.set123.frma.cancer.reorder.st34.MDE.100<-ed.set123.frma.cancer.reorder.st34[which(rownames(ed.set123.frma.cancer.reorder.st34) %in%  TopTableSt34.most.DE.probes.frma[,1]),]

subtype.s41<-c(as.character(check.match[which(check.match[,9]=="st4"),9]),as.character(check.match[which(check.match[,9]=="st1"),9]))
subtype.s34<-c(as.character(check.match[which(check.match[,9]=="st3"),9]),as.character(check.match[which(check.match[,9]=="st4"),9]))

Draw_heatmap(ed.set123.frma.cancer.reorder.st41.MDE.100,"heatmap_allsample_frma_MDE_200_probes_41.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s41)
Draw_heatmap(ed.set123.frma.cancer.reorder.st34.MDE.100,"heatmap_allsample_frma_MDE_200_probes_34.pdf","/media/H_driver/2015/Sophia/Results/",subtype.s34)


Draw_PCA(ed.set123.frma.cancer.reorder,"PCA_58_cancer_sample_frma_all_probes.pdf","/media/H_driver/2015/Sophia/Results/",check.match$subtype)

save(heatmap_wPCA,ed.set123.frma.cancer.reorder,check.match,file="Data_4_heatmap_PCA.RData")

#Use all data to perform GSEA analysis
GenerateFiles4GSEA<-function(ed.set123.frma.cancer,out_dir,out_file_name,subtypeA,subtypeB){

  index.st4<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeA),2])
  index.st1<-which(colnames(ed.set123.frma.cancer) %in% cel.file.name.key.set.symbol.st134[which(cel.file.name.key.set.symbol.st134[,5]==subtypeB),2])

  ed.set123.frma.cancer.st41<-cbind(ed.set123.frma.cancer[,index.st4],ed.set123.frma.cancer[,index.st1])
  cat(subtypeA,"\t",length(index.st4),subtypeB,"\t",length(index.st1),"\t",dim(ed.set123.frma.cancer.st41),"\n")

  ed.set123.frma.cancer.st41.reformat<-cbind(rownames(ed.set123.frma.cancer.st41),rep("NA",dim(ed.set123.frma.cancer.st41)[1]),ed.set123.frma.cancer.st41)
  colnames(ed.set123.frma.cancer.st41.reformat)[c(1,2)]=c("NAME","Description")

  write.table(ed.set123.frma.cancer.st41.reformat,file = paste0(out_dir,out_file_name),quote = FALSE, sep ="\t",row.names = FALSE,col.names = TRUE)

}

GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St41-4_GSEA.xls","st4","st1")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St34-4_GSEA.xls","st3","st4")
GenerateFiles4GSEA(ed.set123.frma.cancer,"/media/H_driver/2015/Sophia/Results/","St13-4_GSEA.xls","st1","st3")

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_New.RData")
savehistory(file="/media/H_driver/Aimin_project/Sophia/Data_set_3_16_2016_all_probes_based_New.Rhistory")
