\name{print.ahazpen}
\alias{print.ahazpen}
\alias{summary.ahazpen}
\title{Print an ahazpen object}
\description{
Print method for fitted penalized semiparametric additive hazards model.
}
\usage{
\method{print}{ahazpen}(x, digits=max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{Fitted \code{ahazpen} object.}
  \item{digits}{Significant digits to print.}
  \item{\dots}{For future methods.}
}
\details{
  The call that produced \code{x} is printed, alongside the
  number of observations, the number of covariates, and details on the
  sequence of penalty parameters.
    }
\seealso{\code{\link{ahazpen}}, \code{\link{predict.ahazpen}}, \code{\link{coef.ahazpen}}.}
\keyword{print}

