\name{predict.tune.ahazpen}
\alias{predict.tune.ahazpen}
\alias{coef.tune.ahazpen}
\title{Prediction methods for tune.ahazpen}
\description{Compute regression coefficient estimates, linear predictor, cumulative hazard function, or
  integrated martingale residuals for a fitted and tuned penalized semiparametric additive
  hazards model.  }
\usage{
\method{predict}{tune.ahazpen}(object, newX,  lambda="lambda.min", \dots)
\method{coef}{tune.ahazpen}(object, \dots)
}
\arguments{
  \item{object}{The result of an \code{ahazpen} fit.}
  \item{newX}{New matrix of covariates at which to do
     predictions. Required for some types of predictions, see \code{\link{predict.ahazpen}}.}
  \item{lambda}{Value of lambda at which predictions are
  to be made. Required for some types of predictions, see
  \code{\link{predict.ahazpen}}. Default is the optimal lambda value
  saved in \code{object}.}
   \item{\dots}{Additional arguments to be passed to
     \code{predict.ahazpen} (usually the type of prediction required).}
}
\details{
  See the details in \code{\link{predict.ahazpen}} for information on
  the available types of predictions.
}


\value{The obejct returned depends on the details in the argument \code{\dots} passed
  to \code{\link{predict.ahazpen}}.}


\seealso{\code{\link{predict.ahazpen}}, \code{\link{ahazpen}}, \code{\link{print.ahazpen}},
  \code{\link{plot.ahazpen}}, \code{\link{predict.ahaz}}, \code{\link{plot.cumahaz}}.}
\examples{
data(sorlie)

set.seed(10101)

# Break ties
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Fit additive hazards regression model w/lasso penalty
cv.fit <- tune.ahazpen(surv, X, dfmax=100, tune="cv")

# Predict coefficients at cv.fit$lambda.min
coef(cv.fit)

# Predict risk score at cv.fit$lambda.min
predict(cv.fit,newX=X,type="lp")

}
\keyword{methods}

 
