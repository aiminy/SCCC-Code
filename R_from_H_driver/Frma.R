biocLite("frmaExampleData")
biocLite("hgu133afrmavecs")
library(frmaExampleData)
library(hgu133afrmavecs)

data(AffyBatchExample)
object <- frma(AffyBatchExample)
