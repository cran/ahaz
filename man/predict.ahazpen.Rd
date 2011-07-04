\name{predict.ahazpen}
\alias{predict.ahazpen}
\alias{coef.ahazpen}
\title{Prediction methods for ahazpen}
\description{Compute regression coefficient estimates, linear predictor, cumulative hazard function, or
  integrated martingale residuals for a fitted penalized semiparametric additive
  hazards model.  }
\usage{
\method{predict}{ahazpen}(object, newX, type=c("coef","lp","residuals","cumhaz"),
        lambda=NULL, \dots)
\method{coef}{ahazpen}(object, \dots)
}
\arguments{
  \item{object}{The result of an \code{ahazpen} fit.}
  \item{newX}{New matrix of covariates at which to do
     predictions. \cr Required unless \code{type="coef"}. }
  \item{lambda}{Value of lambda for at which predictions are
  to be made. This argument is required for \code{type="residuals"} and
  \code{type="cumhaz"}. Since predictions rely on interpolations between lambda values, it is recommended not to use a lambda-value smaller than the minimum of \code{object$lambda}.}
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
   associated with the nonzero penalized regression coefficients for a
   regularization parameter equal to \eqn{lambda}.

   For \code{type="cumhaz"}, an object with S3 class \code{"cumahaz"}
   based on the regression coefficients estimated for a
   regularization parameter equal to \eqn{lambda}, the object containing:
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

# Fit additive hazards regression model w/lasso penalty
fit <- ahazpen(surv, X, dfmax=100)

# Coefficients
beta <- predict(fit,X,lambda=0.08,type="coef")
barplot(as.numeric(beta))

# Linear predictions
linpred <- predict(fit,X,lambda=0.1,type="lp")
riskgrp <- factor(linpred < median(linpred))
plot(survfit(surv~riskgrp))

# Residuals
resid <- predict(fit, X, lambda=0.1, type = "residuals")
par(mfrow = c(1,2))
hist(resid[,1],main=colnames(resid)[1])
hist(resid[,2],main=colnames(resid)[2])

# Cumulative hazard
cumhaz <- predict(fit,X,lambda=0.1,type="cumhaz")
plot(cumhaz)


}
\keyword{methods}

 
