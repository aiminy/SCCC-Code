#load the affy library
biocLite("mogene20sttranscriptcluster.db")
library(gdata)
library(affy)
library(affyio)
library(annotate)
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)

biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

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

cel.file.sample.infor<-as.data.frame(rbind(
cbind(colnames(cancer.data.set123.normalized.st1),rep("st1",length(colnames(cancer.data.set123.normalized.st1)))),
cbind(colnames(cancer.data.set123.normalized.st2),rep("st2",length(colnames(cancer.data.set123.normalized.st2)))),
cbind(colnames(cancer.data.set123.normalized.st3),rep("st3",length(colnames(cancer.data.set123.normalized.st3)))),
cbind(colnames(cancer.data.set123.normalized.st4),rep("st4",length(colnames(cancer.data.set123.normalized.st4))))
))
colnames(cel.file.sample.infor)=c("filename","subtype")

f <- factor(cel.file.sample.infor$subtype)
f
design <- model.matrix(~0+f)
#Differential gene expression analysis for st3 vs st4
probe<-rownames(cancer.data.set123.normalized.st1)
geneSymbol<-getSYMBOL(probe, "hugene10sttranscriptcluster.db")

#Differential gene expression analysis for st2 vs st1


#Differential gene expression analysis for st4 vs st1

save.image(file="/media/H_driver/Aimin_project/Sophia/Data_set_3.RData")
