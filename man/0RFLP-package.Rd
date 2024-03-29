\name{RFLPtools-package}
\alias{RFLPtools-package}
\alias{RFLPtools}
\docType{package}
\title{
Tools To Analyse RFLP-Data
}
\description{
  \pkg{RFLPtools} provides functions to analyse DNA fragment samples 
  (i.e. derived from RFLP-analysis) and standalone BLAST report 
  files (i.e. DNA sequence analysis).
}

\details{
\tabular{ll}{
Package: \tab RFLPtools\cr
Version: \tab 2.0 \cr
Date: \tab 2022-02-07 \cr
Depends: \tab R(>= 4.0.0)\cr
Imports: \tab stats, utils, graphics, grDevices, RColorBrewer\cr
Suggests: \tab knitr, rmarkdown, lattice, MKomics\cr
License: \tab LGPL-3\cr
}}


\author{
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@uni-tuebingen.de},\cr
  Mohammed Aslam Imtiaz,\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de}\cr

\emph{Maintainer:} Matthias Kohl \email{Matthias.Kohl@stamats.de}
}

\references{
  Local Blast download: \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download}
  
  Blast News: \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastNews}
  
  Ian A. Dickie, Peter G. Avis, David J. McLaughlin, Peter B. Reich. 
  Good-Enough RFLP Matcher (GERM) program.
  Mycorrhiza 2003, 13:171-172.

  Flessa, F., Kehl, A., Kohl, M.
  Analysing diversity and community structures using PCR-RFLP: a new software application.
  Molecular Ecology Resources 2013 Jul; 13(4):726-33. 

  Matsumoto, Masaru; Furuya, Naruto; Takanami, Yoichi; Matsuyama, Nobuaki.  
  RFLP analysis of the PCR-amplified 28S rDNA in Rhizoctonia solani. 
  Mycoscience 1996 37:351-356.

  Persoh, D., Melcher, M., Flessa, F., Rambold, G.: 
  First fungal community analyses of endophytic ascomycetes associated with Viscum album ssp. austriacum 
  and itshost Pinus sylvestris. Fungal Biology 2010 Jul;114(7):585-96.

  Poussier, Stephane; Trigalet-Demery, Danielle; Vandewalle, Peggy; Goffinet, Bruno; Luisetti, Jacques; Trigalet, Andre.
  Genetic diversity of Ralstonia solanacearum as assessed by PCR-RFLP of the hrp gene region, AFLP and 16S rRNA sequence 
  analysis, and identification of an African subdivision.
  Microbiology 2000 146:1679-1692.
  
  T. A. Saari,  S. K. Saari, C. D. Campbell, I. J Alexander, I. C. Anderson.
  FragMatch - a program for the analysis of DNA fragment data.
  Mycorrhiza 2007, 17:133-136
}
\examples{
data(RFLPdata)
res <- RFLPdist(RFLPdata)
plot(hclust(res[[1]]), main = "Euclidean distance")

par(mfrow = c(1,2))
plot(hclust(RFLPdist(RFLPdata, nrBands = 3)), cex = 0.7)
RFLPplot(RFLPdata, nrBands = 3, mar.bottom = 6, cex.axis = 0.8)

data(RFLPref)
RFLPrefplot(RFLPdata, RFLPref, nrBands = 6, cex.axis = 0.8)


library(MKomics)
data(BLASTdata)
res <- simMatrix(BLASTdata, sequence.range = TRUE, Min = 500)
myCol <- colorRampPalette(brewer.pal(8, "RdYlGn"))(128)
simPlot(res, col = myCol, minVal = 0, 
        labels = colnames(res), title = "(Dis-)Similarity Plot")
}
\keyword{package}
\keyword{multivariate}
\keyword{cluster}
