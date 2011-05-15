\name{plot.ahaz}
\alias{plot.ahaz}
\title{Plot an ahaz object}
\description{Plot method for a fitted semiparametric additive hazards model;
  plots the Breslow estimate of underlying cumulative hazard function.}
\usage{
\method{plot}{ahaz}(x, \dots)
}
\arguments{
  \item{x}{The result of an \code{ahaz} fit.}
  \item{\dots}{Additional graphical arguments passed to the \code{plot} function.}
}
\details{Calling \code{plot.ahaz} is equivalent to first calling \code{ahaz}, then
  calling \code{predict} with \code{type="cumhaz"}, and finally
  calling \code{plot}.}

\seealso{\code{\link{ahaz}}, \code{\link{predict.ahaz}}, \code{\link{plot.cumahaz}}}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,15:24])

# Fit additive hazards model
fit <- ahaz(surv, X)
plot(fit)
}
\keyword{hplot}