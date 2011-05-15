\name{plot.cumahaz}
\alias{plot.cumahaz}
\title{Plot a cumahaz object}
\description{Plots the Breslow estimate of cumulative hazard function,
  as obtained from  the \code{predict.ahaz}}
\usage{
\method{plot}{cumahaz}(x, \dots)  }

\arguments{
  \item{x}{Result of a call to the \code{predict.ahaz} function with option \code{type="cumhaz"}.}
  \item{\dots}{Additional graphical arguments passed to the \code{plot} function.}
}

\seealso{\code{\link{predict.ahaz}}, \code{\link{predict.ahazpen}}}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,15:24])

# Fit additive hazards regression model
fit <- ahaz(surv, X)

# Cumulative hazard
cumhaz <- predict(fit, type="cumhaz")
plot(cumhaz)
}
\keyword{hplot}
