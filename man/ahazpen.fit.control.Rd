\name{ahazpen.fit.control}
\alias{ahazpen.fit.control}
\title{Controls for ahazpen fitting algorithm}
\description{Controls the numerical algorithm  for fitting
  the penalized semiparametric additive hazards model. This is typically
only used in a call to \code{ahazpen}.}
\usage{
ahazpen.fit.control(thresh=1e-5, maxit=100000, \dots)
}
\arguments{
  \item{thresh}{Declare convergence when the maximal relative change from the
    last iteration is less than \code{thresh}. Default is
    \code{thresh=1e-5}.}
  \item{maxit}{Maximal number passes by the algorithm over the data for
    all values of the regularization parameter lambda. Default is \code{maxit=100000}. }
   \item{\dots}{For future methods.}
}

\value{A list with elements named as the arguments.}

\seealso{\code{\link{ahazpen}}}
\keyword{misc}

 
