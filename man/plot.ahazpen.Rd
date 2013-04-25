\name{plot.ahazpen}
\alias{plot.ahazpen}
\title{Plot an ahazpen object}
\description{Plots regularization paths for fitted penalized semiparametric additive hazards model.}
\usage{
\method{plot}{ahazpen}(x, xvar=c("norm","lambda"), labels=FALSE, df=TRUE,
                ylab="Regression coefficients", xlab=xname,...)}

\arguments{
  \item{x}{The result of an \code{ahazpen} fit.}
  \item{xvar}{Scaling for first axis. Options are the \eqn{L^1}{L1} norm of the vector of
  regression coefficients ("\code{norm}") or the penalty parameter
  on a log scale ("\code{lambda}").}
  \item{labels}{Try to display indices for the regression coefficients in the right-hand
    margin. Default is \code{labels=FALSE}.}
  \item{df}{Display number of nonzero parameters in top margin. Default
    is \code{df=TRUE}.}
  \item{ylab}{Label for y-axis.}
  \item{xlab}{Label for x-axis. The default is either "L1 norm" or
  \eqn{\lambda}{lambda}, depending on \code{xvar}.}
  \item{\dots}{Additional graphical arguments passed to the \code{plot} function.}
}
\seealso{\code{\link{ahazpen}}, \code{\link{print.ahazpen}}, \code{\link{predict.ahazpen}}, \code{\link{coef.ahazpen}}.}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Fit additive hazards regression model
fit <- ahazpen(surv, X, dfmax=50)
par(mfrow=c(1,2)); plot(fit); plot(fit,xvar="lambda")

# With labels only
plot(fit,labels=TRUE,df=FALSE)
}
\keyword{hplot}

 
