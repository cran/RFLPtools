\name{RFLPdist2}
\alias{RFLPdist2}
\title{
  Compute distances for RFLP data.
}
\description{
  If gel image quality is low, faint bands may be disregarded and 
  may lead to wrong conclusions. This function computes the distance 
  between the molecular weights of RFLP samples, including samples 
  containing one or more additional bands. Thus, failures during 
  band detection could be identified. Visualisation of band patterns 
  using this method can be done by \code{\link{RFLPplot}} using the 
  argument \code{nrMissing}.
}
\usage{
RFLPdist2(x, distfun = dist, nrBands, nrMissing, LOD = 0,
          diag = FALSE, upper = FALSE)
}
\arguments{
  \item{x}{ data.frame with RFLP data; see \code{\link{RFLPdata}}.}
  \item{distfun}{ function computing the distance with default \code{dist};
  cf. \code{\link[stats]{dist}}.}
  \item{nrBands}{ samples with number of bands equal to \code{nrBands} are to be considered.}
  \item{nrMissing}{ number of bands that might be missing.}
  \item{LOD}{ threshold for low-bp bands. }
  \item{diag}{ see \code{\link[stats]{dist}} }
  \item{upper}{ see \code{\link[stats]{dist}} }
}
\details{
  For a given number of bands the given distance between the molecular weights
  is computed. It is assumed that a number of bands might be missing. Hence
  all samples with number of bands in nrBands, nrBands+1, ..., nrBands+nrMissing
  are compared. 

  If \code{LOD > 0} is specified, it is assumed that missing bands can only occur for
  molecular weights smaller than \code{LOD}. As a consequence only samples which 
  have \code{nrBands} bands with molecular weight larger or equal to \code{LOD} are
  selected.
  
  For computing the distance between the molecular weight of a sample S1 with x bands 
  and a Sample S2 with x+y bands the distances between the molecular weight of sample S1
  and the molecular weight of all possible subsets of S2 with x bands are computed.
  The distance between S1 and S2 is then defined as the minimum of all these distances.
  
  If \code{LOD > 0} is specified, only all combinations of values below \code{LOD} are 
  considered.
  
  This option may be useful, if gel image quality is low, and the detection of bands 
  is doubtful.
}
\value{
  An object of class \code{"dist"} returned; cf. \code{\link[stats]{dist}}.
}
\references{
  Flessa, F., Kehl, A., Kohl, M.
  Analysing diversity and community structures using PCR-RFLP: a new software application.
  Molecular Ecology Resources 2013 Jul; 13(4):726-33. 
  
  Ian A. Dickie, Peter G. Avis, David J. McLaughlin, Peter B. Reich. 
  Good-Enough RFLP Matcher (GERM) program.
  Mycorrhiza 2003, 13:171-172.
}
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPdata}}, \code{\link{nrBands}}, \code{\link{RFLPdist}}, \code{\link[stats]{dist}} }
\examples{
## Euclidean distance
data(RFLPdata)
nrBands(RFLPdata)
res0 <- RFLPdist(RFLPdata, nrBands = 4)
res1 <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 1)
res2 <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 2)
res3 <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 3)

## assume missing bands only below LOD
res1.lod <- RFLPdist2(RFLPdata, nrBands = 4, nrMissing = 1, LOD = 60)

## hierarchical clustering
par(mfrow = c(2,2))
plot(hclust(res0), main = "0 bands missing")
plot(hclust(res1), main = "1 band missing")
plot(hclust(res2), main = "2 bands missing")
plot(hclust(res3), main = "3 bands missing")

## missing bands only below LOD
par(mfrow = c(1,2))
plot(hclust(res0), main = "0 bands missing")
plot(hclust(res1.lod), main = "1 band missing below LOD")

## Similarity matrix
library(MKomics)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
ord <- order.dendrogram(as.dendrogram(hclust(res1)))
temp <- as.matrix(res1)
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), title = "(Dis-)Similarity Plot")

## missing bands only below LOD
ord <- order.dendrogram(as.dendrogram(hclust(res1.lod)))
temp <- as.matrix(res1.lod)
simPlot(temp[ord,ord], col = rev(myCol), minVal = 0, 
        labels = colnames(temp), title = "(Dis-)Similarity Plot\n1 band missing below LOD")


## or
library(lattice)
levelplot(temp[ord,ord], col.regions = rev(myCol),
          at = do.breaks(c(0, max(temp)), 128),
          xlab = "", ylab = "",
          ## Rotate label of x axis
          scales = list(x = list(rot = 90)),
          main = "(Dis-)Similarity Plot")


## Other distances
res11 <- RFLPdist2(RFLPdata, distfun = function(x) dist(x, method = "manhattan"),
                 nrBands = 4, nrMissing = 1)
res12 <- RFLPdist2(RFLPdata, distfun = corDist, nrBands = 4, nrMissing = 1)
res13 <- RFLPdist2(RFLPdata, distfun = corDist, nrBands = 4, nrMissing = 1, LOD = 60)
par(mfrow = c(2,2))
plot(hclust(res1), main = "Euclidean distance\n1 band missing")
plot(hclust(res11), main = "Manhattan distance\n1 band missing")
plot(hclust(res12), main = "Pearson correlation distance\n1 band missing")
plot(hclust(res13), main = "Pearson correlation distance\n1 band missing below LOD")
}
\keyword{multivariate}
\keyword{cluster}
