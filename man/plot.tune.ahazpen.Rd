\name{plot.tune.ahazpen}
\alias{plot.tune.ahazpen}
\title{Plot a tune.ahazpen object}
\description{Plot, as a function of the penalty parameter, the curve of tuning scores produced when tuning a
  penalized semiparametric additive hazards model.}
\usage{
\method{plot}{tune.ahazpen}(x, df = TRUE, \dots)  }

\arguments{
  \item{x}{The result of a call to \code{tune.ahazpen}.}
     \item{df}{Display number of nonzero parameters in top margin. Default
    is \code{df=TRUE}.}
    \item{\dots}{Additional graphical arguments passed to the \code{plot} function.}
  }

  \details{A plot is produced displaying the tuning score for each value
    of penalty parameter 
    (alongside upper and lower standard deviation curves, if cross-validation
    has been used). The value of
    \code{lambda} which minimizes the estimated tuning score
    is indicated with a dashed vertical line.}
\seealso{\code{\link{ahazpen}}, \code{\link{tune.ahazpen}}, \code{\link{print.tune.ahazpen}}.}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Do 10 fold cross-validation
set.seed(10101)
tune.fit <- tune.ahazpen(surv, X, penalty="lasso",
              dfmax=50, tune = cv.control(nfolds=10))
plot(tune.fit)

}
\keyword{hplot}

 
