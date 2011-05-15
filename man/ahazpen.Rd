\name{ahazpen}
\alias{ahazpen}
\title{Fit penalized semiparametric additive hazards model}
\description{ Fit a  semiparametric additive hazards model via penalized
  estimating equations using, for example, the lasso penalty. The complete regularization path is computed  at a grid
    of values for the penalty parameter lambda. 
  }
\usage{
ahazpen(surv, X, weights,  standardize=TRUE,  penalty=lasso.control(),
        dfmax=nobs-1, lambda, lambda.minf=ifelse(nobs<nvars,0.01,1e-4),
        nlambda=100, penalty.wgt=NULL, keep=NULL, control=list())
}
\arguments{
  \item{surv}{Response in the form of a survival object, as returned by the
  function \code{Surv()} in the package \pkg{survival}. Right-censored
  and counting process format (left-truncation) is supported. Tied
  survival times are not supported.}
  \item{X}{Design matrix. Missing values are not supported.}
  \item{weights}{Optional vector of observation weights. Default is 1 for each observation.}
  \item{standardize}{Logical flag for variable standardization, prior to
    model fitting. Estimates are always returned on
    the original scale. Default is \code{standardize=TRUE}.}
 \item{penalty}{A description of the penalty function to be used for
   model fitting. This can be a character string naming a penalty
   function (currently \code{"lasso"} or stepwise SCAD, \code{"sscad"}) or a call to
   the desired penalty function. \cr See \code{\link{ahazpen.pen.control}} for the available
   penalty functions and advanced options; see also the examples. }
  \item{dfmax}{Limit the maximum number of variables in the
      model. Unless a complete
  regularization path is needed, it is highly
  recommended to initially choose a relatively smaller value of
  \code{dfmax} to substantially reduce computation time. }
    \item{lambda}{An optional user supplied sequence of penalty parameters. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.minf}. }
  \item{nlambda}{The number of \code{lambda} values. Default is \code{nlambda=100}.}
  \item{lambda.minf}{Smallest value of \code{lambda}, as a fraction of
    \code{lambda.max}, the (data-derived) smallest
  value of \code{lambda} for which all regression coefficients are zero. The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs >= nvars}, the default is \code{0.0001},
  close to zero.  When \code{nobs < nvars}, the default is \code{0.01}.} 
\item{penalty.wgt}{A vector of penalty weights for each
  regression coefficient. This is a number that multiplies \code{lambda} to allow
    differential penalization. Can be 0 for some variables, which implies
    no penalization so that the variable is always included in the
    model. Default is 1 for all variables.}
  \item{keep}{A vector of indices of variables which should always be included in
    the model (no penalization). Equivalent to specifying a \code{penalty.wgt} of 0.}
    \item{control}{A list of parameters for controlling the
      model fitting algorithm. The list is passed to \code{\link{ahazpen.fit.control}}.}
  }

\details{

Fits the sequence of models implied by the penalty function
\code{penalty}, the sequence of penalty parameters \code{lambda}, and
  \code{penalty}. 

   For data sets with a very large number of covariates, it is highly
  recommended to only calculate partial paths by specifying a smallish
  value of \code{dmax}.

  The sequence \code{lambda} is computed automatically but can also be
  specified manually. It may happen that the \code{lambda} sequence is
  truncated before \code{nlambda} values of \code{lambda} have been
  used, if \code{dfmax} or a saturated fit is reached. Note that the
  fitting algorithm may exhibit slow convergence, for example when the
  model is close to a saturated fit. In such cases, a warning is issued,
  and one should either increase the value of \code{maxit} via the
  argument \code{control}; or decrease \code{dfmax}.  }

\value{
  An object with S3 class \code{"ahazpen"}.
  \item{call}{The call that produced this object}
  \item{beta}{An \code{
      length(lambda) x nvars} matrix of penalized regression coefficients.}
  \item{lambda}{The sequence of actual \code{lambda} values used.}
  \item{df}{The number of nonzero coefficients for each value of
    \code{lambda}.}
  \item{nobs}{Number of observations.}
  \item{nvars}{Number of covariates.}
  \item{surv}{A copy of the argument \code{survival}.}
   \item{iter}{Maximal number of iterations used during coordinate
     descent loop.}
   \item{penalty.wgt}{The actually used \code{penalty.wgt}.}
   \item{penalty}{An object of class \code{ahaz.pen.control}, as
  specified by \code{penalty}.}
}
\references{
Gorst-Rasmussen, A. & Scheike, T. H. (2011). \emph{Independent screening for
  single-index hazard rate models with ultra-high dimensional features.}
  Technical report, Department of Mathematical Sciences, Aalborg University.

  Leng, C. & Ma, S. (2007). \emph{Path consistent model selection in
 additive risk model via Lasso}. Statistics in Medicine; \bold{26}:3753-3770.

  Martinussen, T. & Scheike, T. H. (2008). \emph{Covariate selection
  for the semiparametric additive risk model.} Scandinavian Journal of
Statistics; \bold{36}:602-619.

 Zou, H. & Li, R. (2008). \emph{One-step sparse estimates in nonconcave
 penalized likelihood models}, Annals of Statistics; \bold{36}:1509-1533.
}
\seealso{\code{\link{print.ahazpen}}, \code{\link{predict.ahazpen}},
 \code{\link{coef.ahazpen}},  \code{\link{plot.ahazpen}},
 \code{\link{tune.ahazpen}}.}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Fit additive hazards regression model
fit1 <- ahazpen(surv, X,penalty="lasso", dfmax=30)
plot(fit1)

# User-specified lambda sequence
lambda <- exp(seq(log(0.30), log(0.1), length = 100))
fit2 <- ahazpen(surv, X, penalty="lasso", lambda = lambda)
plot(fit2)

# Advanced usage - specify details of the penalty function
fit3 <- ahazpen(surv, X,penalty=sscad.control(nsteps=5), dfmax=30)
fit3
fit4 <- ahazpen(surv, X,penalty=lasso.control(alpha=1e-3), dfmax=30)
plot(fit4)

}

\keyword{models}
\keyword{regression}
\keyword{survival}