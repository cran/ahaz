\name{ahaz.partial}
\alias{ahaz.partial}
\title{Partial calculation of estimating quantities used by ahaz}
\description{Partial calculation of the quantities used in the estimating equations for
  ahaz.}
\usage{
ahaz.partial(surv, X, weights, idx)
}

\arguments{
  \item{surv}{Response in the form of a survival object, as returned by the
  function \code{Surv()} in the package \pkg{survival}. Right-censored
  and counting process format (left-truncation) is supported. Tied
  survival times are not supported.}
  \item{X}{Design matrix. Missing values are not supported.}
  \item{weights}{Optional vector of observation weights. Default is 1
    for each observation.}
    \item{idx}{Vector of indices of covariates to use in the calculations.}
}
\details{
  The function is intended mainly for \strong{programming use} when a
  very large number of covariates are considered and direct application of
  \code{ahaz} is unfeasible.
    
  The estimating equations for the semiparametric additive hazards model
  are of the form \eqn{D\beta=d}{D*beta=d} with \eqn{D} a quadratic matrix with
  number of columns equal to the number of covariates. The present
  function returns \code{d[idx]}, \code{D[idx,]}, and \code{B[idx,]};
  the latter a matrix such that \eqn{D^{-1} B
  D^{-1}} estimates the
  covariance matrix of the regression coefficients.
  }
\value{
  A list containing the following elements:
  \item{call}{The call that produced this object.}
  \item{idx}{A copy of the argument \code{idx}.}
  \item{nobs}{Number of observations.}
  \item{nvars}{Number of covariates.}
  \item{d}{Vector of length \code{length(idx)}.}
  \item{D}{Matrix of size \code{length(idx) x nvars}.}
  \item{B}{Matrix of size \code{length(idx) x nvars}.}
}
\seealso{\code{\link{ahaz}}, \code{\link{ahaz.adjust}}.}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Get D for the first 10 covariates only
a<-ahaz.partial(surv,X,idx=1:10)
pD1 <- a$D

# Equivalent to the (slower) solution
b <- ahaz(surv,X)
pD2 <- b$D[1:10,]
max(abs(pD1-pD2))

}
\keyword{methods}

