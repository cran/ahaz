\name{ahazpen}
\alias{ahazpen}
\title{Fit penalized semiparametric additive hazards model}
\description{ Fit a  semiparametric additive hazards model via penalized
  estimating equations using, for example, the lasso penalty. The complete regularization path is computed  at a grid
    of values for the penalty parameter lambda via the method of cyclic
    coordinate descent. 
  }
\usage{
ahazpen(surv, X, weights,  standardize=TRUE,  penalty=lasso.control(),
        nlambda=100, dfmax=nvars, pmax=min(nvars, 2*dfmax),
        lambda.minf=ifelse(nobs < nvars,0.05, 1e-4), lambda,
        penalty.wgt=NULL, keep=NULL, control=list())
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
   \item{nlambda}{The number of \code{lambda} values. Default is
     \code{nlambda=100}.}
   
  \item{dfmax}{Limit the maximum number of variables in the
      model. Unless a complete
  regularization path is needed, it is highly
  recommended to initially choose a relatively smaller value of
  \code{dfmax} to substantially reduce computation time. }
  \item{pmax}{Limit the maximum number of variables to ever be
      considered by the coordinate descent algorithm. }
  \item{lambda.minf}{Smallest value of \code{lambda}, as a fraction of
    \code{lambda.max}, the (data-derived) smallest
  value of \code{lambda} for which all regression coefficients are zero. The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs >= nvars}, the default is \code{0.0001},
  close to zero.  When \code{nobs < nvars}, the default is \code{0.05}.}
    \item{lambda}{An optional user supplied sequence of penalty parameters. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.minf}. A user-specified
      lambda sequence overrides \code{dfmax} but not \code{pmax}.}
\item{penalty.wgt}{A vector of nonnegative penalty weights for each
  regression coefficient. This is a number that multiplies \code{lambda} to allow
    differential penalization. Can be 0 for some variables, which implies
    no penalization so that the variable is always included in the
    model; or \code{Inf} which implies that the variable is never
      included in the model. Default is 1 for all variables.}
  \item{keep}{A vector of indices of variables which should always be included in
    the model (no penalization). Equivalent to specifying a \code{penalty.wgt} of 0.}
    \item{control}{A list of parameters for controlling the
      model fitting algorithm. The list is passed to \code{\link{ahazpen.fit.control}}.}
  }

\details{

Fits the sequence of models implied by the penalty function
\code{penalty}, the sequence of penalty parameters \code{lambda} by
using the very efficient method of cyclic coordinate descent.

For data sets with a very large number of covariates, it is recommended
to only calculate partial paths by specifying a smallish value of
\code{dmax}.

The sequence \code{lambda} is computed automatically by the algorithm
but can also be set (semi)manually by specifying \code{nlambda} or
\code{lambda}. The stability and efficiency of the algorithm is highly
dependent on the grid \code{lambda} values being reasonably dense, and
\code{lambda} (and \code{nlambda}) should be specified accordingly. In
particular, it is not recommended to specify a single or a few lambda
values. Instead, a partial regularization path should be calculated and
the functions \code{\link{predict.ahazpen}} or
\code{\link{coef.ahazpen}} should be used to extract coefficient
estimates at specific lambda values.
}
\value{
  An object with S3 class \code{"ahazpen"}.
  \item{call}{The call that produced this object}
  \item{beta}{An \code{
      nvars x length(lambda)} matrix (in sparse column format, class \code{dgCMatrix}) of penalized regression coefficients.}
  \item{lambda}{The sequence of actual \code{lambda} values used.}
  \item{df}{The number of nonzero coefficients for each value of
    \code{lambda}.}
  \item{nobs}{Number of observations.}
  \item{nvars}{Number of covariates.}
  \item{surv}{A copy of the argument \code{survival}.}
   \item{npasses}{Total number of passes by the fitting algorithm over the data,
  for all lambda values.}
   \item{penalty.wgt}{The actually used \code{penalty.wgt}.}
   \item{penalty}{An object of class \code{ahaz.pen.control}, as
     specified by \code{penalty}.}
      \item{dfmax}{A copy of \code{dfmax}.}
         \item{penalty}{A copy of \code{pmax}.}
}
\references{
  Gorst-Rasmussen A., Scheike T. H. (2012). \emph{Coordinate Descent Methods
  for the Penalized Semiparametric Additive Hazards Model}. Journal of
  Statistical Software, \bold{47}(9):1-17. \url{https://www.jstatsoft.org/v47/i09/}
	
Gorst-Rasmussen, A. & Scheike, T. H. (2011). \emph{Independent screening for
  single-index hazard rate models with ultra-high dimensional features.}
  Technical report R-2011-06, Department of Mathematical Sciences, Aalborg University.

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
fit1
plot(fit1)

# Extend the grid to contain exactly 100 lambda values
lrange <- range(fit1$lambda)
fit2 <- ahazpen(surv, X,penalty="lasso", lambda.minf=lrange[1]/lrange[2])
plot(fit2)

# User-specified lambda sequence
lambda <- exp(seq(log(0.30), log(0.1), length = 100))
fit2 <- ahazpen(surv, X, penalty="lasso", lambda = lambda)
plot(fit2)

# Advanced usage - specify details of the penalty function
fit4 <- ahazpen(surv, X,penalty=sscad.control(nsteps=2))
fit4
fit5 <- ahazpen(surv, X,penalty=lasso.control(alpha=0.1))
plot(fit5)
}

\keyword{models}
\keyword{regression}
\keyword{survival}
