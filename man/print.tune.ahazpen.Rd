\name{print.tune.ahazpen}
\alias{print.tune.ahazpen}
\title{Print a tune.ahazpen object}
\description{Print method for \code{tune.ahazpen} objects.}
\usage{
\method{print}{tune.ahazpen}(x, digits=max(3, getOption("digits") - 3), ...) }

\arguments{
  \item{x}{The result of a call to \code{tune.ahazpen}.}
  \item{digits}{Significant digits in printout.}
  \item{\dots}{Additional print arguments.}
}

\details{
  The call that produced \code{x} is printed, alongside the number of
  penalty parameters used, the value of the
  optimal penalty and the number of non-zero regression coefficients at
  the optimal penalty parameter.
    }
\seealso{\code{\link{ahazpen}}, \code{\link{tune.ahazpen}}, \code{\link{plot.tune.ahazpen}}.}
\keyword{print}
 
