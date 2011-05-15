\name{predict.ahazpen}
\alias{predict.ahazpen}
\alias{coef.ahazpen}
\title{Prediction methods for ahazpen}
\description{Compute regression coefficient estimates, linear predictor, cumulative hazard function, or
  integrated martingale residuals for a fitted penalized semiparametric additive
  hazards model.  }
\usage{
\method{predict}{ahazpen}(object, newX,
    type=c("coef","lp","residuals","cumhaz"), index=NULL, \dots)
\method{coef}{ahazpen}(object, \dots)
}
\arguments{
  \item{object}{The result of an \code{ahazpen} fit.}
  \item{newX}{New matrix of covariates at which to do
     predictions. \cr Required unless \code{type="coef"}. }
  \item{index}{Index of the penalty parameter in \code{object$lambda} at which predictions are
  to be made. Required for \code{type="residuals"} and
  \code{type="cumhaz"}.}
 \item{type}{The type of prediction. Options are the regression coefficients
    ("\code{coef}"), the linear
    predictors ("\code{lp}"),  the (integrated) martingale residuals
    ("\code{residuals}"), or the cumulative hazard  ("\code{cumhaz}")}
   \item{\dots}{For future methods.}
}
\details{
  See the details in \code{\link{predict.ahaz}} for information on
  the different types of predictions.
}


\value{ For \code{type="coef"} and \code{type="lp"}, a
   matrix of regression coefficients, respectively linear predictions for
   each value of the penalty parameter.
   
   For \code{type="residuals"}, a matrix of (integrated) martingale residuals
   associated with the nonzero penalized regression coefficients when the penalty
   parameter in  \code{object$lambda} has index \code{index}.

   For \code{type="cumhaz"}, an object with S3 class \code{"cumahaz"}
   based on the regression coefficients estimated for the penalty
   parameter in \code{object$lambda} with index \code{index}, and containing:
   \item{time}{Jump times for the cumulative hazard estimate.}
   \item{cumhaz}{The cumulative hazard estimate.}
   \item{event}{Status at jump times (1 corresponds to death, 0
     corresponds to entry/exit).}

}


\seealso{\code{\link{ahazpen}}, \code{\link{print.ahazpen}},
  \code{\link{plot.ahazpen}}, \code{\link{predict.ahaz}}, \code{\link{plot.cumahaz}}.}
\examples{
data(sorlie)

set.seed(10101)

# Break ties
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Fit additive hazards regression model
fit <- ahazpen(surv, X)

# Fit additive hazards regression model with lasso penalty
fit <- ahazpen(surv, X,penalty="lasso")

# Coefficients
beta <- predict(fit,X,type="coef")
barplot(beta[50,])

# Linear predictions
linpred <- predict(fit,X,type="lp")
riskgrp <- factor(linpred[15,] < median(linpred[15,]))
plot(survfit(surv~riskgrp))

# Residuals
resid <- predict(fit, X, index=15, type = "residuals")
par(mfrow = c(1,2))
hist(resid[,1],main=colnames(resid)[1])
hist(resid[,2],main=colnames(resid)[2])

# Cumulative hazard
cumhaz <- predict(fit,X,index=15,type="cumhaz")
plot(cumhaz)


}
\keyword{methods}

 
