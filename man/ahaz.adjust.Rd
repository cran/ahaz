\name{ahaz.adjust}
\alias{ahaz.adjust}
\title{Adjusted univariate association measures from ahaz}
\description{Fast calculation of univariate association measures  in the
  semiparametric additive risk model, adjusted for user-specified covariates}
\usage{
ahaz.adjust(surv, X, weights, idx, method=c("coef", "z", "crit"))
}

\arguments{
  \item{surv}{Response in the form of a survival object, as returned by the
  function \code{Surv()} in the package \pkg{survival}. Right-censored
  and counting process format (left-truncation) is supported. Tied
  survival times are not supported.}
  \item{X}{Design matrix. Missing values are not supported.}
  \item{weights}{Optional vector of observation weights. Default is 1
    for each observation.}
    \item{idx}{Vector specifying the indices of the covariates to adjust
      for.}
    \item{method}{The type of adjusted association measure to calculate. See details.}
}
\details{
   The function is intended mainly for \strong{programming use} and
  screening purposes, when a very large number of covariates are considered and direct application of
  \code{ahaz} is unfeasible.
  
  Running this function is equivalent to running \code{ahaz} with
  design matrix \code{cbind(X[,i],X[,idx])} for each column \code{X[,i]} in
  \code{X}. By utilizing basic matrix identities, \code{ahaz.adjust}
  runs many times faster.
  
  The following univariate association measures are currently implemented:
  \itemize{
    \item  \code{method="z"},  \eqn{Z}-statistics, obtained from a
  fitted \code{ahaz} model.
    \item  \code{method="coef"}, regression coefficients, obtained from a
  fitted \code{ahaz} model.
    \item \code{method="crit"}, the increase in
    the natural loss function of the
    semiparametric additive hazards model when the covariate is included
  in the model.
  }
}
\value{
  A list containing the following elements:
  \item{call}{The call that produced this object.}
  \item{idx}{A copy of the argument \code{idx}.}
  \item{adj}{Adjusted association statistic, as specified by \code{method}. Entries with index in
  \code{idx} are set to \code{NA}.}
}
\seealso{\code{\link{ahaz}}, \code{\link{ahaz.partial}}, \code{\link{ahazisis}}.}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Adjust for first 10 covariates
idx <- 1:10
a <- ahaz.adjust(surv,X,idx=idx)

# Compare with (slower) solution
b <- apply(X[,-idx],2,function(x){coef(ahaz(surv,cbind(x,X[,idx])))[1]})
plot(b,a$adj[-idx])
}
\keyword{models}
\keyword{regression}
\keyword{survival}

 
