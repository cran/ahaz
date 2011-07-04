\name{predict.ahaz}
\alias{predict.ahaz}
\alias{coef.ahaz}
\alias{vcov.ahaz}
\alias{residuals.ahaz}
\title{Prediction methods for ahaz}
\description{Compute regression coefficients, linear predictor, cumulative hazard function, or
  integrated martingale residuals for a fitted semiparametric additive hazards model. }
\usage{
\method{predict}{ahaz}(object, newX, type=c("coef", "lp",
       "residuals", "cumhaz"), beta=NULL, \dots)
\method{coef}{ahaz}(object, \dots)
\method{vcov}{ahaz}(object, \dots)
\method{residuals}{ahaz}(object, \dots)
}
\arguments{
  \item{object}{The result of an \code{ahaz} fit.}
   \item{newX}{Optional new matrix of covariates at which to do
     predictions. Currently only supported for \code{type="lp"}.}
  \item{type}{Type of prediction. Options are the regression coefficients
    ("\code{coef}"), the linear
    predictor ("\code{lp}"),  the martingale residuals
    ("\code{residuals}"), or the cumulative hazard
    ("\code{cumhaz}"). See the details.}
  \item{beta}{Optional vector of regression coefficients. If unspecified,
    the regression coefficients derived from \code{object} are used.}
  \item{\dots}{For future methods.}
  }
  \details{
    The Breslow estimator of the baseline cumulative hazard is described
    in Lin & Ying (1994).

    The regression coefficients \eqn{\beta_0}{beta_0} in the semiparametric additive hazards
    model are obtained as the
    solution \eqn{\hat{\beta}}{beta} to a quadratic  system of linear equations
    \eqn{D\beta=d}{D*beta = d}. The
    (integrated) martingale residuals \eqn{\epsilon_i}{epsilon_i} for
    \eqn{i=1,...,n} are vectors, of length
    corresponding to the number of covariates, so that
    \deqn{D(\hat{\beta}-\beta_0) -d \approx \epsilon_1+\cdots+\epsilon_n}{D*(beta-beta_0) - d
      ~ epsilon_1 + ... + epsilon_n}
    The residuals estimate integrated
martingales and are
    asymptotically distributed as mean-zero IID multivariate Gaussian. They can be used to derive a sandwich-type variance
    estimator for regression coefficients (implemented in
    \code{summary.ahaz} when \code{robust=TRUE} is specified).  They can moreover be used for implementing consistent standard error
    estimation under clustering; or for implementing resampling-based
    inferential methods.

    See Martinussen & Scheike (2006), Chapter 5.4 for details.  }
    \value{

    For \code{type="coef"} and \code{type="lp"}, a vector of
    predictions.

    For \code{type="coef"}, a matrix of (integrated) martingale
    residuals, with number of columns corresponding to the number of
    covariates.

   For \code{type="cumhaz"}, an object with S3 class \code{"cumahaz"}
 consisting of:
 \item{time}{Jump times for the cumulative hazard estimate.}
 \item{cumhaz}{The cumulative hazard estimate.}
  \item{event}{Status at jump times (1 corresponds to death, 0
 corresponds to entry/exit).}
}
\references{Martinussen, T. & Scheike, T. H. &  (2006). \emph{Dynamic Regression Models for Survival Data.} Springer.}
\seealso{\code{\link{ahaz}},  \code{\link{summary.ahaz}}, \code{\link{plot.cumahaz}}.}
\examples{
data(sorlie)

set.seed(10101)

# Break ties
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,15:24])

# Fit additive hazards regression model
fit <- ahaz(surv, X)

# Parameter estimates
coef(fit)

# Linear predictor, equivalent to X\%*\%coef(fit)
predict(fit,type="lp")

# Cumulative baseline hazard
cumahaz <- predict(fit, type="cumhaz")

# Residuals - model fit
resid <- predict(fit, type = "residuals")
# Decorrelate, standardize, and check QQ-plots
stdres <- apply(princomp(resid)$scores,2,function(x){x/sd(x)})
par(mfrow = c(2,2))
for(i in 1:4){
  qqnorm(stdres[,i])
  abline(c(0,1))
}

# Residuals - alternative variance estimation
resid <- residuals(fit)
cov1 <- summary(fit)$coef[,2]
invD <- solve(fit$D)
Best<-t(resid)\%*\%resid
cov2 <- invD \%*\% Best \%*\% invD
# Compare with (nonrobust) SEs from 'summary.ahaz'
plot(cov1, sqrt(diag(cov2)),xlab="Nonrobust",ylab="Robust")
abline(c(0,1))
}
\keyword{methods}

