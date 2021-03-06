\name{nrBands}
\alias{nrBands}
\title{
  Function to compute number of bands.
}
\description{
  Computes groups based on the number of bands per sample in a RFLP data set.
  Each group comprises RFLP-samples with equal number of bands.
}
\usage{
nrBands(x)
}
\arguments{
  \item{x}{ data.frame with RFLP data; see \code{\link{RFLPdata}}.}
}
\details{
  The function computes groups based on the number of bands per sample 
  in a RFLP data set. Each group comprises RFLP-samples with equal 
  number of bands.
}
\value{ Number of bands per RFLP-samples. }
\references{
  Flessa, F., Kehl, A., Kohl, M.
  Analysing diversity and community structures using PCR-RFLP: a new software application.
  Molecular Ecology Resources 2013 Jul; 13(4):726-33. 
}
\author{ 
  Fabienne Flessa \email{Fabienne.Flessa@uni-bayreuth.de},\cr
  Alexandra Kehl \email{Alexandra.Kehl@uni-tuebingen.de},\cr
  Matthias Kohl \email{Matthias.Kohl@stamats.de} 
}
%\note{}
\seealso{ \code{\link{RFLPdata}}, \code{\link{RFLPdist2}}, \code{\link[stats]{dist}} }
\examples{
data(RFLPdata)
nrBands(RFLPdata)
}
\keyword{multivariate}
