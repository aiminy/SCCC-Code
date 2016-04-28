#load the affy library
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

# Read in the CEL files in the directory, then normalize the data
setwd("/media/H_driver/2015/Sophia/FTE-HGSC Gene Expression/Microarray_Set1_FTE_HGSC_LCM")
data.set1 <- ReadAffy()
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

colnames(data.set1.normalized)
colnames(data.set2.normalized)
colnames(data.set3.normalized)

length(colnames(data.set1.normalized))
length(colnames(data.set2.normalized))
length(colnames(data.set3.normalized))

unique(colnames(data.set1.normalized))
length(unique(colnames(data.set1.normalized)))
length(unique(colnames(data.set2.normalized)))
length(unique(colnames(data.set3.normalized)))

data.set1.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/Microarray_Set1_FTE_HGSC_LCM/CASE_DESCRIPTION_SET1.xlsx")
data.set2.sample.info<-read.xls("/home/aiminyan/Downloads/GSE28044_George.xlsx")
data.ER.PR.sample.info<-read.xls("/media/H_driver/2015/Sophia/FTE-HGSC\ Gene\ Expression/ER_PR_Subgroup_PT2_v2.xlsx")
head(data.ER.PR.sample.info)

head(data.set1.sample.info)
data.set1.sample.info[,2]
data.set1.normalize
data.set2.normalized

grep(2367,unique(data.set1.sample.info[,2]))

unique(data.set2.sample.info[,2])
head(data.ER.PR.sample.info)

dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==2),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==3),])
dim(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==4),])

MapSample2CELfile<-function(data.ER.PR.sample.info,data.set1.normalized,sample_type){
  
cat(dim(data.ER.PR.sample.info),dim(data.set1.normalized),sample_type,"\n")
  
num.sample1=length(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==1),][,1])
cel_file_index_4_sample<-array()

for(i in 1:num.sample1) {
tmp<-grep(data.ER.PR.sample.info[which(data.ER.PR.sample.info[,4]==sample_type),][i,1],colnames(data.set1.normalized))
cel_file_index_4_sample<-c(cel_file_index_4_sample,tmp)
}
cel_file_index_4_sample<-cel_file_index_4_sample[-1]

colnames(data.set1.normalized)[cel_file_index_4_sample]

re<-data.set1.normalized[,cel_file_index_4_sample]

return(re)

}

#ER-PR-
ER.PR.sample1.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,1)
ER.PR.sample1.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,1)
ER.PR.sample1.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,1)

#ER-PR+
ER.PR.sample2.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,2)
ER.PR.sample2.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,2)
ER.PR.sample2.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,2)

#ER+PR-
ER.PR.sample3.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,3)
ER.PR.sample3.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,3)
ER.PR.sample3.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,3)

#ER+PR+
ER.PR.sample4.from.set1<-MapSample2CELfile(data.ER.PR.sample.info,data.set1.normalized,4)
ER.PR.sample4.from.set2<-MapSample2CELfile(data.ER.PR.sample.info,data.set2.normalized,4)
ER.PR.sample4.from.set3<-MapSample2CELfile(data.ER.PR.sample.info,data.set3.normalized,4)

head(ER.PR.sample1.from.set1)
head(ER.PR.sample1.from.set2)
head(ER.PR.sample1.from.set3)

head(ER.PR.sample3.from.set1)
head(ER.PR.sample3.from.set2)
head(ER.PR.sample3.from.set3)

head(data.set1.sample.info)
head(data.set2.sample.info)
head(data.ER.PR.sample.info)

data.set1.sample.info[,2]

data.set2.sample.info[,8:9:10]

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
dim(cancer.data.set123.normalized)

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

dim(cancer.data.set123.normalized.st1)
dim(cancer.data.set123.normalized.st2)
dim(cancer.data.set123.normalized.st3)
dim(cancer.data.set123.normalized.st4)

head(cancer.data.set123.normalized.st1)
head(cancer.data.set123.normalized.st2)
head(cancer.data.set123.normalized.st3)
head(cancer.data.set123.normalized.st4)

colnames(cancer.data.set123.normalized.st1)
colnames(cancer.data.set123.normalized.st2)
colnames(cancer.data.set123.normalized.st3)
colnames(cancer.data.set123.normalized.st4)

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

colnames(cel.file.sample.infor.no.2)=c("filename","subtype")
f.st134 <- factor(cel.file.sample.infor.no.2$subtype)
design.st134 <- model.matrix(~0+f.st134)
colnames(design.st134) <- levels(f.st134)

cancer.data.st134<-cbind(cancer.data.set123.normalized.st1,
cancer.data.set123.normalized.st3,
cancer.data.set123.normalized.st4)
head(cancer.data.st134)
GeneSym.all <- getSYMBOL(rownames(cancer.data.st134), "hgu133plus2.db")

ndata<-cancer.data.st134
geneSymbol=GeneSym.all
tempdata.byGSym = data.frame(ndata, Symbol = geneSymbol)

tempdata.byGSym.2<-tempdata.byGSym[-which(is.na(tempdata.byGSym$Symbol)),]

rownames(tempdata.byGSym.2) = NULL
data.byGSym = ddply(tempdata.byGSym.2, c("Symbol"),function(h)
  summary = apply(h,2,max)
)

#class(data.byGSym)

#which(is.na(data.byGSym$Symbol))
rownames(data.byGSym)=data.byGSym$Symbol

data.byGSym.2<-data.byGSym[,-61]

dim(data.byGSym.2)

#Use all probes

fit.st134 <- lmFit(cancer.data.st134, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.all<-topTable(fit2.st134,coef=1,n=54674)

TopTableSt34.100<-topTable(fit2.st134,coef=1,adjust="fdr",n=100)

genenames <- as.character(rownames(TopTableSt34.54675))

length(genenames)

genenames

annotation(ed.normalized.set1)
annotation(ed.normalized.set2)
annotation(ed.normalized.set3)
library("hgu133plus2.db")

map <- getpAnnMa("ENTREZID", "hgu133plus2", load=TRUE, type=c("env", "db"))
class(map)
ll <- getEG(genenames,"hgu133plus2.db")
GeneSym <- getSYMBOL(genenames, "hgu133plus2.db")

Probe.gene.sym<-(cbind(genenames,GeneSym))

Probe.gene.sym

length(which(is.na(Probe.gene.sym[,2])))


tab <- data.frame(GeneSym, TopTableSt34.100)
tab <- data.frame(rownames(tab), tab) 

colnames(tab)[1] <- c("Probe ID") 
ll <- list(ll)
htmlpage(ll, filename="/media/H_driver/2015/Sophia/St34-4.html", title="HTML report", 
         othernames=tab, table.head=c("Locus ID",colnames(tab)), table.center=TRUE, digits=6)

colnames(fit2.st134)
topTable(fit2.st134,coef=1)
topTable(fit2.st134,coef=1,adjust="fdr")

topTable(fit2.st134,coef=2)
topTable(fit2.st134,coef=2,adjust="fdr")

results.st134 <- decideTests(fit2.st134,p.value=0.99)

summary(results.st134)
vennDiagram(results.st134)

fit2.st134$genes

unigeneTopTableSt34 <- topTable(fit2.st134,coef=1,n=20,genelist=genelist)
unigeneTopTableSt41 <- topTable(fit2.st134,coef=2,n=20,genelist=genelist)

library(xtable)
xtableUnigeneSt34 <- xtable(unigeneTopTableSt34,display=c("s","s","s","s","g","g","g","e","e","g","g"))
xtableUnigeneSt41 <- xtable(unigeneTopTableSt41,display=c("s","s","s","s","g","g","g","e","e","g","g"))

cat(file="/media/H_driver/2015/Sophia/St34.html","<html>\n<body>")
print.xtable(xtableUnigeneSt34,type="html",file="/media/H_driver/2015/Sophia/St34.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St34.html","</body>\n</html>",append=TRUE)

cat(file="/media/H_driver/2015/Sophia/St41.html","<html>\n<body>")
print.xtable(xtableUnigeneSt41,type="html",file="/media/H_driver/2015/Sophia/St41.html",append=TRUE)
cat(file="/media/H_driver/2015/Sophia/St41.html","</body>\n</html>",append=TRUE)


#Use gene name
temp = apply(data.byGSym.2,2,as.numeric)
rownames(temp) = rownames(data.byGSym.2)
data.byGSym.2 = temp


fit.st134.gene <- lmFit(data.byGSym.2, design.st134)

cont.matrix.st134 <- makeContrasts(st34="st3-st4",st41="st4-st1",levels=design.st134)
fit2.st134  <- contrasts.fit(fit.st134.gene, cont.matrix.st134)
fit2.st134  <- eBayes(fit2.st134)

str(fit2.st134)

TopTableSt34.gene<-topTable(fit2.st134,coef=1,n=54674)
dim(TopTableSt34.gene)
head(TopTableSt34.gene)
hist(TopTableSt34.gene[,4])

TopTableSt41.gene<-topTable(fit2.st134,coef=2,n=54674)
dim(TopTableSt41.gene)
head(TopTableSt41.gene)
hist(TopTableSt41.gene[,4])


heatmap_wPCA = function(Data, output_heatmap, output_pca, out_dir, g_level = NULL){
  hmcol<-rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
  if(is.null(g_level)){
    type_level = 1:ncol(Data)
    col_level = "black"
  }else{
    type_level = 16
    TEMP = factor(g_level)
    uniq_label =  levels(TEMP)
    levels(TEMP) = hmcol[ceiling(seq(length.out=length(levels(TEMP)),from=1,to=256))]
    col_level = as.character(TEMP)
    uniq_col = levels(TEMP)
  }
  
  Data.hc = hclust(dist(Data), method="average")
  rowInd <- order.dendrogram(as.dendrogram(Data.hc))
  pdf(file=paste0(out_dir,"/",output_heatmap), width=10, height=12)
  heatmap.2(Data[rowInd,], col=hmcol, Rowv = F, dendrogram = "column", scale="row",labRow = NA,key=TRUE, keysize=0.55, symkey=FALSE,density.info="none", trace="none",cexCol=1,margins=c(8,12))
  dev.off()
  
  pdf(file=paste0(out_dir,"/",output_pca), width=12, height=12)
  Data.pca = prcomp(t(Data))
  with(data.frame(Data.pca$x), scatter3D(PC1, PC2, PC3, colvar = NULL, type="h",
                                         ticktype = "detailed", bty="b2", cex=1,
                                         xlab="PC 1",	ylab="PC 2",zlab="PC 3", theta = 40, phi = 40, pch=type_level,
                                         col=col_level,
                                         main = "Principal component analysis")
  )
  legend("topright", legend = uniq_label, pch=type_level, 
         col = uniq_col, 
         cex=1, inset=c(0.02))
  with(data.frame(Data.pca$x), text3D(x=PC1, y=PC2, 
                                      z=PC3, colnames(Data), col = "black", add=TRUE, colkey = FALSE, cex=0.5)
  )
  dev.off()  
}

SampleType = factor(gsub("(FTY).[1-3]","\\1",colnames(data.byGSym)))

SampleType<-cel.file.sample.infor.no.2$subtype

heatmap_wPCA(data.byGSym.2, "heatmap_allsample.pdf","PCA_allsample.pdf","/media/H_driver/2015/Sophia/", SampleType)

# clustering plot

#fit#Differential gene expression analysis for st3 vs st4
#probe<-rownames(cancer.data.set123.normalized.st1)
#geneSymbol<-getSYMBOL(probe, "hugene10sttranscriptcluster.db")
#Differential gene expression analysis for st2 vs st1
#Differential gene expression analysis for st4 vs st1
save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")