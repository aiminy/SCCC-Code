model.file<-dir(pattern="bam_mm_model.r")

#dir.name="/media/H_driver/2016/Morey_project/Peak_chip_seq/Shift75_summits/"
#file.name=dir("/media/H_driver/2016/Morey_project/Peak_chip_seq/Shift75_summits/",pattern="summits")
  
dir.name="/media/H_driver/2016/Morey_project/Peak_chip_rm_control/"
file.name=dir("/media/H_driver/2016/Morey_project/Peak_chip_rm_control/",pattern="summits")

peak.files<-paste0(dir.name,file.name)
name.sample<-gsub("__bam_mm_shift_75_2_summits.bed","",gsub(dir.name,"",peak.files))


peak.files.list<-as.list(peak.files)
names(peak.files.list)=name.sample


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peak.files.list, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList[[1]], xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")
tagHeatmap(tagMatrixList[[1]], xlim=c(-3000, 3000), title="7",color="red")

names(peak.files.list[[1]])

for(i in 1:18)
{
source(model.file[i])
}

source("http://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
biocLite("clusterProfiler")

require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)

files <- getSampleFiles()
print(files)
class(files)
names(files)

peak <- readPeakFile(files[[4]])
covplot(peak, weightCol="V5")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")



peakAnno <- annotatePeak(peak.files[1], tssRegion=c(-3000, 3000), TxDb=txdb)
peakAnno.data<-as.data.frame(peakAnno)

dim(peakAnno.data)

head(peakAnno.data)

class(peakAnno)
csAnno(peakAnno)

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)


upsetplot(peakAnno, vennpie=TRUE)

plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")

as.data.frame(peakAnno)


getGEOspecies()



promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peak.files.list, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")

tagHeatmap(tagMatrixList[[7]], xlim=c(-3000, 3000), title="7",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="8",color="red")
tagHeatmap(tagMatrixList[[13]], xlim=c(-3000, 3000), title="13",color="red")
tagHeatmap(tagMatrixList[[12]], xlim=c(-3000, 3000), title="12",color="red")




tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")

tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[5]], xlim=c(-1000, 1000), title="H3K4me3",color="red")
tagHeatmap(tagMatrixList[[8]], xlim=c(-3000, 3000), title="H3K4me3",color="red")






