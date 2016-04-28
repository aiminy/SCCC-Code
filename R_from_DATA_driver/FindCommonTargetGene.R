Gene.H3K4me3<-read.csv("H3K4me3vsInp.tss.R1.gene_name.csv")
Gene.H3K27me3<-read.csv("H3K27me3vsInp.tss.R1.gene_name.csv")
Gene.H3K4me3.H3K27me3<-read.csv("H3K4me3--H3K27me3vsInp.tss.R1.gene_name.csv")
Gene.H3K27me3.H3K4me3<-H3K4me3<-read.csv("K3K27me3--H3K4me3vsInp.tss.R1.gene_name.csv")

dim(Gene.H3K4me3)
dim(Gene.H3K27me3)

head(Gene.H3K4me3[order(Gene.H3K4me3[,1]),])
head(Gene.H3K27me3[order(Gene.H3K27me3[,1]),])

gene.common<-intersect(Gene.H3K4me3[,1],Gene.H3K27me3[,1])

gene.common.re.chip<-intersect(gene.common,Gene.H3K4me3.H3K27me3[,1])

length(gene.common.re.chip)

head(Gene.H3K4me3)

head(Gene.H3K27me3)

dir()



H3K27me3.gene<-read.table("H3K27me3-2.txt",header = F)
H3K4me3.gene<-read.table("H3K4me3-2.txt",header = F)
H3K4me3.H3K27me3.gene<-read.table("H3K4me3-H3K27me3-2.txt")
H3K27me3.H3K4me3.gene<-read.table("H3K27me3-H3K4me3-2.txt")



gene.common<-intersect(H3K27me3.gene[,1],H3K4me3.gene[,1])

length(gene.common.2)

gene.common.1<-intersect(gene.common,H3K4me3.H3K27me3.gene[,1])

gene.common.2<-intersect(gene.common,H3K27me3.H3K4me3.gene[,1])

