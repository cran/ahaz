\name{summary.ahaz}
\alias{summary.ahaz}
\title{Summarize an ahaz object}
\description{Produces a summary of a fitted semiparametric additive
  hazards model.}
\usage{
\method{summary}{ahaz}(object, ...)  }
\arguments{
  \item{object}{The result of an \code{ahaz} fit.}
  \item{...}{For future methods.}
}
\value{
  An object with S3 class \code{"summary.ahaz"}.
  \item{call}{The call that produced this object.}
  \item{coefficients}{Vector of regression coefficients.}
  \item{cov}{Estimated covariance matrix of regression coefficients.}
  \item{nobs}{Number of observations.}
  \item{nvars}{Number of covariates}
  \item{waldtest}{Vector of quantities from a Wald test.}
  \item{univar}{Logical: summarizing univariate
  regressions (option \code{univariate} in \code{ahaz})?}
}

\seealso{\code{\link{ahaz}}, \code{\link{plot.ahaz}}}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,15:25])

# Fit additive hazards model
fit1 <- ahaz(surv, X)
summary(fit1)
}
\keyword{methods}

 
