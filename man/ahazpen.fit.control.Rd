\name{ahazpen.fit.control}
\alias{ahazpen.fit.control}
\title{Controls for ahazpen fitting algorithm}
\description{Controls the numerical algorithm  for fitting
  the penalized semiparametric additive hazards model. This is typically
only used in a call to \code{ahazpen}.}
\usage{
ahazpen.fit.control(thresh=1e-4, maxit=2500, \dots)
}
\arguments{
  \item{thresh}{Declare convergence when the maximal relative change from the
    last iteration is less than \code{thresh}. Default is \code{thresh=1e-4}.}
  \item{maxit}{Maximal number of coordinate descent iterations before the algorithm exits
    with a warning. Default is \code{maxit=2500}. If very
    large models are considered or the penalized model is closed to
    saturated, it may be necessary to increase \code{maxit}.}
   \item{\dots}{For future methods.}
}

\value{A list with elements named as the arguments.}

\seealso{\code{\link{ahazpen}}}
\keyword{misc}

 
