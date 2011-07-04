\name{ahaz}
\alias{ahaz}
\title{Fit semiparametric additive hazards model}
\description{Fit a semiparametric additive hazards regression
  model. Right-censored and left-truncated survival data are supported.}
\usage{ahaz(surv, X, weights, univariate=FALSE, robust=FALSE)  }

\arguments{
  \item{surv}{Response in the form of a survival object, as returned by the
  function \code{Surv()} in the package \pkg{survival}. Right-censoring
  and left-truncation is supported. Tied survival times are not supported.}
  \item{X}{Design matrix. Missing values are not supported.}
  \item{weights}{Optional vector of observation weights. Default is 1 for each observation.}
  \item{univariate}{Fit all univariate models instead of the joint
    model. Default is \code{univar =
      FALSE}.}
    \item{robust}{Robust calculation of variance. Default is \code{robust = FALSE}.
  }
}
\details{
  The semiparametric additive hazards model specifies a hazard function
  of the form:
  \deqn{h(t) = h_0(t) + \beta' Z_i}{h(t) = h_0(t) + beta' * Z_i}
  for \eqn{i=1,\ldots,n}{i=1,...,n} where \eqn{Z_i} is the vector of covariates,
  \eqn{\beta}{beta} the vector of regression coefficients and \eqn{h_0} is an
  unspecified baseline hazard. The semiparametric additive hazards model
  can be viewed as an
  additive analogue of the well-known Cox proportional hazards
  regression model.

  Estimation is based on the estimating equations of Lin & Ying (1994).

  The option \code{univariate} is intended for screening purposes in
  data sets with a large number of covariates. It is substantially faster than the
  standard approach of combining \code{ahaz} with
  \code{apply}, see the examples.
  }


\value{
  An object with S3 class \code{"ahaz"}.
  \item{call}{The call that produced this object.}
  \item{nobs}{Number of observations.}
  \item{nvars}{Number of covariates.}
  \item{D}{A \code{nvars x nvars}  matrix (or vector of length
  \code{nvars} if \code{univar =
      TRUE}).}
  \item{d}{A vector of length \code{nvars}; the
  regression coefficients equal \code{solve(D,d)}.}
  \item{B}{An \code{nvars x nvars} matrix such that \eqn{D^{-1} B
  D^{-1}} estimates the
  covariance matrix of the regression coefficients. If \code{robust=FALSE}
  then \code{B} is estimated using an asymptotic approximation; if
  \code{robust=TRUE} then \code{B} is estimated from residuals, see \code{\link[=residuals.ahaz]{residuals}}.}
\item{univariate}{Is \code{univariate=TRUE}?}
\item{data}{Formatted version of original data (for internal use).}
  \item{robust}{Is \code{robust=TRUE}?}
}
\references{Lin, D.Y. & Ying, Z. (1994). \emph{Semiparametric analysis of
  the additive risk model.} Biometrika; \bold{81}:61-71.}
\seealso{
\code{\link{summary.ahaz}}, \code{\link{predict.ahaz}},
\code{\link{plot.ahaz}}.
The functions \code{\link[=coef.ahaz]{coef}},
\code{\link[=vcov.ahaz]{vcov}}, \code{\link[=residuals.ahaz]{residuals}}.
}

\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,15:24])

# Fit additive hazards model
fit1 <- ahaz(surv, X)
summary(fit1)

# Univariate models
X <- as.matrix(sorlie[,3:ncol(sorlie)])
fit2 <- ahaz(surv, X, univariate = TRUE)
# Equivalent to the following (slower) solution
beta <- apply(X,2,function(x){coef(ahaz(surv,x))})
plot(beta,coef(fit2))

}

\keyword{models}
\keyword{regression}
\keyword{survival}