\name{tune.ahazpen}
\alias{tune.ahazpen}
\title{Choice of penalty parameter in ahazpen}
\description{Tuning of penalty parameters for  the  penalized semiparametric
  additive hazards model via cross-validation - or via non-stochastic
  procedures, akin to BIC for likelihood-based models.
  }
\usage{
tune.ahazpen(surv, X, weights, standardize=TRUE, penalty=lasso.control(),
             tune=cv.control(), dfmax=nvars, lambda, ...)
}

\arguments{
  \item{surv}{Response in the form of a survival object, as returned by the
  function \code{Surv()} in the package \pkg{survival}. Right-censored
  and counting process format (left-truncation) is supported. Tied
  survival times are not supported.}
  \item{X}{Design matrix. Missing values are not supported.}
  \item{weights}{Optional vector of observation weights. Default is 1 for each observation.}
  \item{standardize}{Logical flag for variable standardization, prior to
    model fitting. Parameter estimates are always returned on
    the original scale. Default is \code{standardize=TRUE}.}
 \item{penalty}{A description of the penalty function to be used for
   model fitting. This can be a character string naming a penalty
   function (currently \code{"lasso"} or stepwise SCAD, \code{"sscad"})
  or it can be a call to the penalty function. Default is
  \code{penalty=lasso.control()}. See \code{\link{ahazpen.pen.control}} for the available
   penalty functions and advanced options; see also the examples.  }
    \item{dfmax}{Limit the maximum number of covariates included in the
      model. Default is \code{nvars=nobs-1}. Unless a complete
  regularization path is needed, it is \strong{highly}
  recommended to initially choose a relatively smaller value of
  \code{dfmax} to reduce computation time and memory usage.}
\item{lambda}{An optional user supplied sequence of penalty parameters. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.min}.}
  \item{tune}{A description of the tuning method to be used. This can be
  a character string naming a tuning control
   function (currently \code{"cv"} or \code{"bic"}) or a call to the
  tuning control function. Default is
  5-fold cross-validation, 
   \code{tune=cv.control()}, see \code{\link{ahaz.tune.control}} for more
  options. See also the examples.
     }
  \item{...}{Additional arguments to be passed to \code{ahazpen}, see
  \code{\link{ahazpen}} for options.}
}
\details{
  The function performs an initial penalized fit based on the
  penalty supplied in \code{penalty} to obtain a  sequence of
  penalty parameters. Subsequently, it selects among these an optimal penalty parameter based on
  the tuning control function described in \code{tune}, see  \code{\link{ahaz.tune.control}}.
}
\value{
  An object with S3 class \code{"tune.ahazpen"}.
  \item{call}{The call that produced this object.}
  \item{lambda}{The actual sequence of \code{lambda} values used.}
  \item{tunem}{The tuning score for each value of \code{lambda} (mean
  cross-validated error, if \code{tune=cv.control()}).}
  \item{tunesd}{Estimate of the cross-validated standard error, if \code{tune=cv.control()}.}
  \item{tunelo}{Lower curve = \code{tunem-tunemsd}, if \code{tune=cv.control()}.}
  \item{tuneup}{Upper curve = \code{tunem+tunemsd}, if \code{tune=cv.control()}.}
  \item{lambda.min}{Value of \code{lambda} for which \code{tunem} is
    minimized.}
  \item{df}{Number of non-zero coefficients at each value of
    \code{lambda}.}
  \item{tune}{The selected \code{tune} of S3 class
    \code{"ahaz.tune.control"}.}
  \item{penalty}{The selected \code{penalty} of S3 class \code{"ahazpen.pen.control"}.}
  \item{foldsused}{Folds actually used, if \code{tune=cv.control()}.}
}
\references{
Gorst-Rasmussen, A. & Scheike, T. H. (2011). \emph{Independent screening for
  single-index hazard rate models with ultra-high dimensional features.}
  Technical report R-2011-06, Department of Mathematical Sciences, Aalborg University.
}
\seealso{\code{\link{ahaz.tune.control}}, \code{\link{plot.tune.ahazpen}}, \code{\link{ahazpen}}.}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Training/test data
set.seed(20202)
train <- sample(1:nrow(sorlie),76)
test <- setdiff(1:nrow(sorlie),train)

# Run cross validation on training data
set.seed(10101)
cv.las <- tune.ahazpen(surv[train,], X[train,],dfmax=30)
plot(cv.las)

# Check fit on the test data
testrisk <- predict(cv.las,X[test,],type="lp")
plot(survfit(surv[test,]~I(testrisk<median(testrisk))),main="Low versus high risk")

# Advanced example, cross-validation of one-step SCAD
# with initial solution derived from univariate models.
# Since init.sol is specified as a function, it is
# automatically cross-validated as well
scadfun<-function(surv,X,weights){coef(ahaz(surv,X,univariate=TRUE))}
set.seed(10101)
cv.ssc<-tune.ahazpen(surv[train,],X[train,],
                     penalty=sscad.control(init.sol=scadfun),
                     tune=cv.control(rep=5),dfmax=30)
# Check fit on test data
testrisk <- predict(cv.ssc,X[test,],type="lp")
plot(survfit(surv[test,]~I(testrisk<median(testrisk))),main="Low versus high risk")
}

\keyword{methods}
\keyword{survival}

 
