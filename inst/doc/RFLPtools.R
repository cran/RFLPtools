## -----------------------------------------------------------------------------
library(RFLPtools)

## -----------------------------------------------------------------------------
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "AZ091016_report.txt")
RFLP1 <- read.rflp(file = filename)
str(RFLP1)

RFLP2 <- RFLPqc(RFLP1, rm.band1 = FALSE) # identical to RFLP1
identical(RFLP1, RFLP2)

RFLP3 <- RFLPqc(RFLP1)
str(RFLP3)

RFLP4 <- RFLPqc(RFLP1, rm.band1 = TRUE, QC.rm = TRUE)
str(RFLP4)

## -----------------------------------------------------------------------------
data(RFLPdata)
res <- RFLPdist(RFLPdata)
names(res) ## number of bands
str(res$"6")

## -----------------------------------------------------------------------------
res1 <- RFLPdist(RFLPdata, distfun = function(x) dist(x, method = "manhattan"))
res2 <- RFLPdist(RFLPdata, distfun = function(x) dist(x, method = "maximum"))
str(res[[1]])
str(res1[[1]])
str(res2[[1]])

## -----------------------------------------------------------------------------
library(MKomics)
res3 <- RFLPdist(RFLPdata, distfun = corDist)
str(res3$"9")

## ---- fig.height=7, fig.width=7-----------------------------------------------
plot(hclust(res[[1]]), main = "Euclidean distance")
plot(hclust(res1[[1]]), main = "Manhattan distance")
plot(hclust(res2[[1]]), main = "Maximum distance")
plot(hclust(res3[[1]]), main = "Pearson correlation distance")

## -----------------------------------------------------------------------------
clust4bd <- hclust(res[[2]])
cgroups50 <- cutree(clust4bd, h=50)
cgroups50

## ---- fig.height=7, fig.width=7-----------------------------------------------
library(RColorBrewer)
library(MKomics)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
ord <- order.dendrogram(as.dendrogram(hclust(res[[1]])))
temp <- as.matrix(res[[1]])
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), title = "(Dis-)Similarity Plot")

## ---- fig.height=7, fig.width=7-----------------------------------------------
library(lattice)
print(levelplot(temp[ord,ord], col.regions = rev(myCol),
          at = do.breaks(c(0, max(temp)), 128),
          xlab = "", ylab = "",
          ## Rotate labels of x-axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot"))

## -----------------------------------------------------------------------------
## Euclidean distance
data(RFLPdata)
data(RFLPref)
nrBands(RFLPdata)
res0 <- RFLPdist(RFLPdata, nrBands = 9)
res1 <- RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 1)
res2 <- RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 2)
res3 <- RFLPdist2(RFLPdata, nrBands = 9, nrMissing = 3)

## ---- fig.height=7, fig.width=7-----------------------------------------------
plot(hclust(res0), main = "0 bands missing")
plot(hclust(res1), main = "1 band missing")
plot(hclust(res2), main = "2 bands missing")
plot(hclust(res3), main = "3 bands missing")

## ---- fig.height=7, fig.width=7-----------------------------------------------
RFLPdata.lod <- RFLPlod(RFLPdata, LOD = 60)
par(mfrow = c(1, 2))
RFLPplot(RFLPdata, nrBands = 4, ylim = c(40, 670))
RFLPplot(RFLPdata.lod, nrBands = 4, ylim = c(40, 670))
title(sub = "After applying RFLPlod")

## ---- fig.height=7, fig.width=7-----------------------------------------------
res0 <- RFLPdist(RFLPdata, nrBands = 4)
res1.lod <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 1, LOD = 60)
ord <- order.dendrogram(as.dendrogram(hclust(res1.lod)))
temp <- as.matrix(res1.lod)
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), 
        title = "(Dis-)Similarity Plot\n1 band missing below LOD")

## ---- fig.height=7, fig.width=7-----------------------------------------------
RFLPrefplot(RFLPdata, RFLPref, nrBands = 9, cex.axis = 0.8)

## ---- eval=FALSE--------------------------------------------------------------
#  system("blastn -query Testquery -db Testdatabase -outfmt 6 -out out.txt")

## ---- eval=FALSE--------------------------------------------------------------
#  ## -outfmt 6
#  test.res <- read.blast(file = "out.txt")

## ---- eval=FALSE--------------------------------------------------------------
#  ## -outfmt 10
#  test.res <- read.blast(file = "out.csv", sep = ",")

## -----------------------------------------------------------------------------
Dir <- system.file("extdata", package = "RFLPtools") # input directory 
filename <- file.path(Dir, "BLASTexample.txt")
BLAST1 <- read.blast(file = filename)
str(BLAST1)

## -----------------------------------------------------------------------------
data(BLASTdata)

## -----------------------------------------------------------------------------
res <- simMatrix(BLASTdata)

## -----------------------------------------------------------------------------
res1 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 100, Max = 450)
res2 <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 500)

## ---- fig.height=7, fig.width=7-----------------------------------------------
library(MKomics)
simPlot(res2, col = myCol, minVal = 0, cex.axis = 0.5,
        labels = colnames(res2), title = "(Dis-)Similarity Plot")

## ---- fig.height=7, fig.width=7-----------------------------------------------
library(lattice)
txt <- trellis.par.get("add.text")
txt$cex <- 0.5
trellis.par.set("add.text" = txt)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
print(levelplot(res2, col.regions = myCol,
          at = do.breaks(c(0, max(res2)), 128),
          xlab = "", ylab = "", 
          ## Rotate labels of x axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot"))

## -----------------------------------------------------------------------------
res.d <- sim2dist(res2)

## ---- fig.height=7, fig.width=7-----------------------------------------------
## hierarchical clustering
plot(hclust(res.d), cex = 0.7)

## -----------------------------------------------------------------------------
sessionInfo()

