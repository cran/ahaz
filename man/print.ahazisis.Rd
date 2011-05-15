\name{print.ahazisis}
\alias{print.ahazisis}
\title{Print an ahazisis object}
\description{
Print method for sure independence screening based on the
semiparametric additive hazards model.
}
\usage{
\method{print}{ahazisis}(x, digits=max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{Fitted \code{ahazisis} object.}
  \item{digits}{Significant digits to print.}
  \item{\dots}{For future methods.}
}
\details{
  The call that produced \code{x} is printed, alongside the
  number of covariates initially recruited, the number of covariates
  finally recruited (if applicable) and the number of iterations (if applicable).
    }
\seealso{\code{\link{ahazisis}}}
\keyword{print}
