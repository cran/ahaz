\name{print.summary.ahaz}
\alias{print.summary.ahaz}
\title{Print a summary.ahaz object}
\description{Produces a printed summary of a fitted semiparametric additive
  hazards model.}
\usage{
\method{print}{summary.ahaz}(x, digits=max(getOption("digits") - 3, 3),
      signif.stars=getOption("show.signif.stars"), \dots)  }

\arguments{
  \item{x}{The result of a call to \code{summary.ahaz}.}
  \item{digits}{Significant digits to print.}
  \item{signif.stars}{Show stars to highlight small p-values.}
  \item{\dots}{For future methods.}
}
\seealso{\code{\link{summary.ahaz}}, \code{\link{ahaz}}, \code{\link{plot.ahaz}}.}
\keyword{print}

 
