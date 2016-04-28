source("http://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
biocLite("clusterProfiler")
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(clusterProfiler)

#load data base
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#read file
peak.files<-paste0("/media/H_driver/2016/Morey_project/Peak_chip_seq/Shift75_summits/",dir("./Shift75_summits"))
name.sample<-gsub("__bam_mm_shift_75_2_summits.bed","",gsub("/media/H_driver/2016/Morey_project/Peak_chip_seq/Shift75_summits/","",peak.files))
peak.files.list<-as.list(peak.files)
names(peak.files.list)=name.sample

peak.1 <- readPeakFile(peak.files[[1]])

tagMatrix.1 <- getTagMatrix(peak.1, windows=promoter)
plotAvgProf(tagMatrix.1, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peak.files.list, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList[[1]],xlim=c(-3000, 3000),facet = "column")

plotAvgProf(tagMatrixList[[7]], xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")

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






