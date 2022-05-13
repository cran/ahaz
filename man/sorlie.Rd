\name{sorlie}
\alias{sorlie}
\docType{data}
\title{Sorlie gene expressions}
\description{
 Dataset containing 549 gene expression measurement, exit time and exit status in a study of breast cancer
 among 115 women.
}
\usage{data(sorlie)}
\format{
  \describe{
    \item{time}{Time to exit.}
    \item{status}{Status at exit (censoring = 0, event = 1).}
    \item{X1,...,X549}{Gene expression measurements.} 
  }
}
\references{
  Soerlie T., et al. (2003). \emph{Repeated observation
     of breast tumor subtypes in independent gene expression data
     sets}. Proc Natl Acad Sci \bold{100}:8418-8423
}

\examples{
  data(sorlie)
}
\keyword{datasets}
