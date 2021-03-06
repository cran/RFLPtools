\name{RFLPdata}
\alias{RFLPdata}
\docType{data}
\title{ Example data set for RFLP data }
\description{
  This is an example data set for RFLP data. 
}
\usage{data(RFLPdata)}
\format{
  A data frame with 737 observations on the following four variables
  \describe{
    \item{\code{Sample}}{character: sample identifier. }
    \item{\code{Band}}{integer: band number. }
    \item{\code{MW}}{integer: molecular weight. }
    \item{\code{Gel}}{character: gel identifier. }
  }
}
\details{
  The molecular weight was determined using the software package Gene Profiler 4.05 
  (Scanalytics Inc.) for DNA fragment analysis and genotyping, and exported to a 
  text file.
}
\source{
  The data set was generated by F. Flessa.
}
\references{
  Flessa, F., Kehl, A., Kohl, M.
  Analysing diversity and community structures using PCR-RFLP: a new software application.
  Molecular Ecology Resources 2013 Jul; 13(4):726-33. 
}
\examples{
data(RFLPdata)
str(RFLPdata)
}
\keyword{datasets}
