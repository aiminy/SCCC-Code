source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(limma)
library(gplots)
library(VennDiagram)
library(UpSetR)
library(org.Mm.eg.db)
library(biomaRt)

dir.name="/media/H_driver/2016/Morey_project/GeneAnno/"
name.pattern="_from_broadPeak_gene_name.txt"

file.name=dir("/media/H_driver/2016/Morey_project/GeneAnno/",pattern=name.pattern)

peak.files<-paste0(dir.name,file.name)
name.sample<-gsub(name.pattern,"",gsub(dir.name,"",peak.files))

peak.files.list<-as.list(peak.files)
names(peak.files.list)=name.sample

#peak.files.list.2<-peak.files.list[-which(names(peak.files.list)=="Peak2Gene")]

ProcessPeak2Gene<-function(x){
  Name.gene<-read.table(x,header = F)
  y<-Name.gene[,1]
  return(y)
}

peak.2.gene<-lapply(peak.files.list,ProcessPeak2Gene)
names(peak.2.gene)

#venn(peak.2.gene[c(8,5,9,7)])

venn.plot <- venn.diagram(
  x = peak.2.gene[c(8+1,5+1,9+1,7+1)],
  filename = "case1_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue", "green"),
  alpha = 0.50,
  label.col = c(rep("white",15)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue", "green"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)

venn.plot <- venn.diagram(
  x = peak.2.gene[c(3+1,5+1,6+1,4+1)],
  filename = "case2_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue", "green"),
  alpha = 0.50,
  label.col = c(rep("white",15)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue", "green"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)

venn.plot <- venn.diagram(
  x = peak.2.gene[c(14+2,11+2,2+1,13+2)],
  filename = "case3_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue", "green"),
  alpha = 0.50,
  label.col = c(rep("white",15)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue", "green"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)

#venn(peak.2.gene[c(12,14,13)])

venn.plot <- venn.diagram(
  x = peak.2.gene[c(12+2,14+2,13+2)],
  filename = "case4_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue"),
  alpha = 0.50,
  label.col = c(rep("white",7)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)

venn.plot <- venn.diagram(
  x = peak.2.gene[c(10+2,14+2,11+2,15+2)],
  filename = "case5_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue", "green"),
  alpha = 0.50,
  label.col = c(rep("white",15)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue", "green"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)

venn.plot <- venn.diagram(
  x = peak.2.gene[c(10+2,1,12+2,14+2,13+2)],
  filename = "case6_new_call.tiff",
  height = 3000, 
  width = 3500, 
  resolution = 1000,
  col = "black",
  lty = "dotted",
  lwd = 1,
  fill = c("red", "orange", "blue", "green","darkorchid4"),
  alpha = 0.50,
  label.col = c(rep("white",31)),
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  #cat.col = c("red", "orange", "blue", "green","darkorchid4"),
  cat.cex = 0.5,
  #cat.just = list(c(-1, -1,1,1,1)),
  #cat.just=c(1,2,3,4,5),
  #rotation.degree = 30,
  cat.pos = 0,
  cat.dist = 0.05,
  #rotation.degree = 30,
  #rotation = TRUE,
  cat.fontfamily = "serif"
)


venn.plot <- venn.diagram(
  x = peak.2.gene[c(12+2,13+2,14+2)],
  filename = "case7_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue"),
  alpha = 0.50,
  label.col = c(rep("white",7)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)

#Cbx7 vx Cbx7--H3K27me3

venn.plot <- venn.diagram(
  x = peak.2.gene[c(1,2)],
  filename = "/media/DATA/ngsplot/case_Cbx7_vs_Cbx7_H3K27me3_new_call.tiff",
  height = 3000, 
  width = 3500, 
  resolution = 1000,
  col = "black",
  lty = "dotted",
  lwd = 1,
  fill = c("red","blue"),
  alpha = 0.50,
  label.col = c(rep("white",3)),
  cex = 0.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red","blue"),
  cat.cex = 0.5,
  #cat.just = list(c(-1, -1,1,1,1)),
  #cat.just=c(1,2,3,4,5),
  #rotation.degree = 30,
  cat.pos = 0.5,
  cat.dist = 0.05,
  #rotation.degree = 30,
  #rotation = TRUE,
  cat.fontfamily = "serif"
)



#common.targte.gene.H3K4cme3.H3K27me3.2<-toupper(common.targte.gene.H3K4cme3.H3K27me3)
#length(common.targte.gene.H3K4me3.H3K27me3)
#library(biomaRt)
#annn.mm10<-read.csv("/media/H_driver/Annotation/mm10/genes_table_02052016.csv")
#length(which(annn.mm10[,1] %in% common.targte.gene.H3K4cme3.H3K27me3.2))
#write.table(annn.mm10[which(annn.mm10[,1] %in% common.targte.gene.H3K4cme3.H3K27me3.2),3],file="H3K4cme3_H3K27me3_common.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)
peak.2.gene

common.targte.gene.H3K4cme3.H3K27me3<-intersect(peak.2.gene[[8+1]],peak.2.gene[[5+1]])
write.table(common.targte.gene.H3K4cme3.H3K27me3,file="/media/DATA/ngsplot/H3K4cme3_H3K27me3_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.H2AK119ub1.H3K27me3<-intersect(peak.2.gene[[3+1]],peak.2.gene[[5+1]])
write.table(common.targte.gene.H2AK119ub1.H3K27me3,file="/media/DATA/ngsplot/H2AK119ub1_H3K27me3_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.Suz12<-peak.2.gene[[14+2]]
write.table(common.targte.gene.Suz12,file="/media/DATA/ngsplot/Suz12.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.Rybp.Suz12<-intersect(peak.2.gene[[12+2]],peak.2.gene[[14+2]])
write.table(common.targte.gene.Rybp.Suz12,file="/media/DATA/ngsplot/Rybp_Suz12_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.Ring1B.Suz12<-intersect(peak.2.gene[[10+2]],peak.2.gene[[14+2]])
write.table(common.targte.gene.Ring1B.Suz12,file="/media/DATA/ngsplot/Ring1B_Suz12_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)


common.targte.gene.Ring1B.Cbx7<-intersect(peak.2.gene[[10+2]],peak.2.gene[[1]])
common.targte.gene.Ring1B.Cbx7.Rybp<-intersect(common.targte.gene.Ring1B.Cbx7,peak.2.gene[[12+2]])
common.targte.gene.Ring1B.Cbx7.Rybp.Suz12<-intersect(common.targte.gene.Ring1B.Cbx7.Rybp,peak.2.gene[[14+2]])
write.table(common.targte.gene.Ring1B.Cbx7.Rybp.Suz12,file="/media/DATA/ngsplot/Ring1B_Cbx7_Rybp_Suz12_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.Rybp.Rybp.Suz12<-intersect(peak.2.gene[[12+2]],peak.2.gene[[13+2]])
write.table(common.targte.gene.Rybp.Rybp.Suz12,file="/media/DATA/ngsplot/Rybp_Rybp_Suz12_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.Rybp.Rybp.Suz12.Suz12<-intersect(common.targte.gene.Rybp.Rybp.Suz12,peak.2.gene[[14+2]])
Suz12.gene.but.not.in.common.targte.gene.Rybp.Rybp.Suz12.Suz12<-setdiff(peak.2.gene[[14+2]],common.targte.gene.Rybp.Rybp.Suz12)

new.Suz12.gene<-setdiff(common.targte.gene.Rybp.Rybp.Suz12,common.targte.gene.Rybp.Rybp.Suz12.Suz12)
write.table(common.targte.gene.Rybp.Rybp.Suz12.Suz12,file="/media/DATA/ngsplot/Rybp_Rybp_Suz12_Suz12_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)
write.table(new.Suz12.gene,file="/media/DATA/ngsplot/new_Suz12_gene.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

write.table(peak.2.gene[[12+2]],file="/media/DATA/ngsplot/Rybp.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)
write.table(peak.2.gene[[13+2]],file="/media/DATA/ngsplot/Rybp_Suz12.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

common.targte.gene.Cbx7.Cbx7.H3K27me3<-intersect(peak.2.gene[[1]],peak.2.gene[[2]])
#length(common.targte.gene.Cbx7.Cbx7.H3K27me3)
write.table(common.targte.gene.Cbx7.Cbx7.H3K27me3,file="/media/DATA/ngsplot/Cbx7_H3K27me3_common_2.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

gene.in.cbx7.but.not.in.H3K27me3<-setdiff(peak.2.gene[[1]],common.targte.gene.Cbx7.Cbx7.H3K27me3)
write.table(gene.in.cbx7.but.not.in.H3K27me3,file="/media/DATA/ngsplot/Cbx7_only_but_not_in_Cbx7_H3K27me3.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

gene.in.H3K27me3.but.not.in.cbx7<-setdiff(peak.2.gene[[2]],common.targte.gene.Cbx7.Cbx7.H3K27me3)
write.table(gene.in.H3K27me3.but.not.in.cbx7,file="/media/DATA/ngsplot/Cbx7_H3K27me3_only_but_not_in_Cbx7.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

sample.dir="/media/DATA/metaseq/metaseq-Lluis-data/"
sample.bam.2<-dir(sample.dir,pattern = "2016-03-20-1-2_S9_.bam")[1]
sample.bam.18<-dir(sample.dir,pattern = "2016-03-23-2-18_S7_.bam")[1]
input.bam<-dir(sample.dir,pattern = "2016-03-20-1-1_S5_.bam")[1]

cat(file="Cbx7_based_on_cbx7_H3K27me3_common.txt",paste0(sample.dir,sample.bam.2,":",sample.dir,input.bam," Cbx7_H3K27me3_common_2.txt"," \"Cbx7\""))
cat(file="Cbx7_H3K27me3_based_on_cbx7_H3K27me3_common.txt",paste0(sample.dir,sample.bam.18,":",sample.dir,input.bam," Cbx7_H3K27me3_common_2.txt"," \"Cbx7-H3K27me3\""))

cat(file="Cbx7_based_on_cbx7_only.txt",paste0(sample.dir,sample.bam.2,":",sample.dir,input.bam," Cbx7_only_but_not_in_Cbx7_H3K27me3.txt"," \"Cbx7\""))
cat(file="Cbx7_H3K27me3_based_on_cbx7_only.txt",paste0(sample.dir,sample.bam.18,":",sample.dir,input.bam," Cbx7_only_but_not_in_Cbx7_H3K27me3.txt"," \"Cbx7-H3K27me3\""))

cat(file="Cbx7_based_on_H3K27me3_only.txt",paste0(sample.dir,sample.bam.2,":",sample.dir,input.bam," Cbx7_H3K27me3_only_but_not_in_Cbx7.txt"," \"Cbx7\""))
cat(file="Cbx7_H3K27me3_based_on_H3K27me3_only.txt",paste0(sample.dir,sample.bam.18,":",sample.dir,input.bam," Cbx7_H3K27me3_only_but_not_in_Cbx7.txt"," \"Cbx7-H3K27me3\""))




Case.all<-read.table("NewBashRun.sh",header=FALSE)
Case.all.2<-read.table("All_config.txt",header=FALSE)
label.case<-c(rep("case1",4),rep("case2",4),rep("case3",4),rep("case4",3),rep("case5",4),rep("case6",5),rep("case7",3))
label.case.2<-cbind(Case.all,label.case)

dim(Case.all)

label.case.3<-c(rep(as.character(Case.all.2[1,1]),4),rep(as.character(Case.all.2[2,1]),4),rep(as.character(Case.all.2[3,1]),4),
                rep(as.character(Case.all.2[4,1]),3),rep(as.character(Case.all.2[5,1]),4),rep(as.character(Case.all.2[6,1]),5),
                rep(as.character(Case.all.2[10,1]),3))
  
re1<-cbind(label.case.2,label.case.3)

label.case.4<-c(rep(as.character(Case.all.2[1,1]),4),rep(as.character(Case.all.2[2,1]),4),rep(as.character(Case.all.2[3,1]),4),
                rep(as.character(Case.all.2[4,1]),3),rep(as.character(Case.all.2[5,1]),4),rep(as.character(Case.all.2[6,1]),5),
                rep(as.character(Case.all.2[7,1]),2),rep(as.character(Case.all.2[9,1]),1))

re2<-cbind(label.case.2,label.case.4)


OutputFunction<-function(case,re1){
  rre<-re1[which(re1[,3]==case),c(1,4,2)]

#  print(rre)
  
  rre1c<-paste0("/media/DATA/metaseq/",rre[,1],":/media/DATA/metaseq/metaseq-Lluis-data/2016-03-20-1-1_S5_.bam",sep="")
  rre2c<-as.character(rre[,2])
  rre3c<-paste0("\"", rre[,3], "\"",sep="")
  
  rrre<-cbind(rre1c,rre2c,rre3c)
  
#  print(rrre)
  
  write.table(rrre,file=paste0("/media/DATA/ngsplot/config",".",case,".txt"),quote=FALSE,
            row.names = FALSE,col.names = FALSE)
  
}

OutputFunction("case1",re1)
OutputFunction("case2",re1)
OutputFunction("case3",re1)
OutputFunction("case4",re1)
OutputFunction("case5",re1)
OutputFunction("case6",re1)
OutputFunction("case7",re1)

dev.off()

OutputFunction("case1",re2)
OutputFunction("case2",re2)
OutputFunction("case3",re2)
OutputFunction("case4",re2)
OutputFunction("case5",re2)
OutputFunction("case6",re2)
OutputFunction("case7",re2)

#head(annn.mm10)
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
# read in the file
#mart = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl")
#getBM(attributes=c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=common.targte.gene.H3K4me3.H3K27me3, mart=mart)

hm1<-read.table("/media/DATA/ngsplot/heatmap_2_1_5kb_4_18_2016/hm1.txt",header=TRUE)
hm2<-read.table("/media/DATA/ngsplot/heatmap_2_2_5kb_4_18_2016/hm1.txt",header=TRUE)
hm3<-read.table("/media/DATA/ngsplot/heatmap_2_3_5kb_4_18_2016/hm1.txt",header=TRUE)
hm4<-read.table("/media/DATA/ngsplot/heatmap_2_4_5kb_4_18_2016/hm1.txt",header=TRUE)

tail(cbind(as.character(hm1[,5]),as.character(hm2[,5]),as.character(hm3[,5]),as.character(hm4[,5])))

hm1[1:10,1:5]
class(hm1)
dim(hm1)
heatmap(
as.matrix(hm1[,5:101]))

image(t(as.matrix(hm1[,5:101])))
image(t(as.matrix(hm2[,5:101])))
image(t(as.matrix(hm3[,5:101])))
image(t(as.matrix(hm4[,5:101])))
axis(1, at=1, labels=1, lwd=1, lwd.ticks=1)

heatmap.2(as.matrix(hm1[,5:101]), dendrogram = "none", Rowv = FALSE, Colv = FALSE) 
heatmap.2(as.matrix(hm2[,5:101]), dendrogram = "none", Rowv = FALSE, Colv = FALSE)
heatmap.2(as.matrix(hm3[,5:101]), dendrogram = "none", Rowv = FALSE, Colv = FALSE)
heatmap.2(as.matrix(hm4[,5:101]), dendrogram = "none", Rowv = FALSE, Colv = FALSE)

hm3.max<-read.table("/media/DATA/ngsplot/heatmap_max_rm_3_5kb_4_18_2016/hm1.txt",header=TRUE)
hm3.total<-read.table("/media/DATA/ngsplot/heatmap_total_rm_3_5kb_4_18_2016/hm1.txt",header=TRUE)

dim(hm3.max)
dim(hm3.total)

image(t(as.matrix(hm3.max[,5:105][order(apply(hm3.max[,5:105],1,max)),])))

image(t(as.matrix(hm3.max[,5:105][order(apply(hm3.max[,5:105],1,sum)),])))

image2(t(as.matrix(hm3.max[,5:105][which.max(apply(hm3.max[,5:105],1,sum)),])))

range(apply(hm3.max[,5:105],1,max))
range(apply(hm3.max[,5:105],1,sum))


#Process Cbx7 and Cbx7--H3K27me3
suffix.name="_5kb_gene_name.txt"
name.pattern="Cbx7*5kb_gene_name.txt"
file.name=dir("/media/H_driver/2016/Morey_project/GeneAnno/",pattern=glob2rx(name.pattern))

peak.files<-paste0(dir.name,file.name)
name.sample<-gsub(suffix.name,"",gsub(dir.name,"",peak.files))

peak.files.list<-as.list(peak.files)
names(peak.files.list)=name.sample

peak.2.gene.Cbx7.Cbx7.H3K27me3<-lapply(peak.files.list,ProcessPeak2Gene)
names(peak.2.gene.Cbx7.Cbx7.H3K27me3)


venn.plot <- venn.diagram(
  x = peak.2.gene[c(12+2,13+2,14+2)],
  filename = "case7_new_call.tiff",
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red", "orange", "blue"),
  alpha = 0.50,
  label.col = c(rep("white",7)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red", "orange", "blue"),
  cat.cex = 0.8,
  cat.fontfamily = "serif"
)


venn.plot <- venn.diagram(
  x = peak.2.gene.Cbx7.Cbx7.H3K27me3[c(1,2)],
  filename = "/media/DATA/ngsplot/case_Cbx7_vs_Cbx7_H3K27me3_original_call.tiff",
  #height = 3000, 
  #width = 3500, 
  #resolution = 1000,
  col = "black",
  lty = "dotted",
  lwd = 2,
  fill = c("red","blue"),
  alpha = 0.50,
  label.col = c(rep("white",3)),
  cex = 1,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("red","blue"),
  cat.cex = 0.8,
  #cat.just = list(c(-1, -1,1,1,1)),
  #cat.just=c(1,2,3,4,5),
  #rotation.degree = 30,
  #cat.pos = 0.5,
  #cat.dist = 0.05,
  #rotation.degree = 30,
  #rotation = TRUE,
  cat.fontfamily = "serif"
)







savehistory(file="/media/DATA/ngsplot/Chip-Seq.Rhistory")
save.image(file="/media/DATA/ngsplot/Chip-Seq.RData")
