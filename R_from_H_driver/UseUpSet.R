movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
mutations <- read.csv( system.file("extdata", "mutations.csv", package = "UpSetR"), header=T, sep = ",")

dim(mutations)
class(mutations)

upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

upset(mutations[1:10,1:5], sets = c("TTN", "PTEN", "TP53","EGFR"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

names(peak.2.gene)


cbind(peak.2.gene[1],peak.2.gene[2])


z <- list(I=c("A","B","C","D"),II=c("A","B"),III=c("B","C","D"))
z.2<-unique(do.call(c,z))

d.test<-matrix(c(1,1,1,1,1,1,0,0,0,1,1,1),4.3)
rownames(d.test)<-z.2
colnames(d.test)<-names(z)

upset(as.data.frame(d.test), sets.bar.color = "#56B4E9",order.by = "freq", empty.intersections = "on")


peak.2.gene.2<-do.call(rbind,peak.2.gene)

peak.2.gene.2[,1:10]