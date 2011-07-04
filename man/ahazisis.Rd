\name{ahazisis}
\alias{ahazisis}
\title{Independent screening for the semiparametric
  additive hazards model}
\description{Fast and scalable model selection for the semiparametric additive hazards
  model via univariate screening combined with penalized regression.
  }
\usage{
ahazisis(surv, X, weights, standardize=TRUE,
        nsis=floor(nobs/1.5/log(nobs)), do.isis=TRUE,
        maxloop=5, penalty=sscad.control(), tune=cv.control(),
        rank=c("FAST","coef","z","crit"))
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
 \item{nsis}{Number of covariates to recruit initially. If
  \code{do.isis=TRUE}, then this is also the maximal number of variables
  that the algorithm will recruit.  Default is
   \code{nsis=floor(nobs/log(nobs)/1.5)}}. 
   \item{do.isis}{Perform iterated independent screening?}
   \item{maxloop}{Maximal number of iterations of the algorithm if \code{do.isis=TRUE}.}
   \item{rank}{Method to use for (re)recruitment of variables. See details.}
 \item{penalty}{A description of the penalty function to be used for
   the variable selection part. This can be a character string naming a penalty
   function  (currently \code{"lasso"} or stepwise SCAD, \code{"sscad"})
  or a call to the penalty function. Default is
  \code{penalty=sscad.control()}. See \code{\link{ahazpen}} and
  \code{\link{ahazpen.pen.control}} for more options and examples.}
  \item{tune}{A description of the tuning method to be used for the
  variable selection part. This can be
  a character string naming a tuning control
   function (currently \code{"cv"} or \code{"bic"}) or a call to the tuning control function. Default is
   \code{tune=cv.control()}.  See \code{\link{ahaz.tune.control}} for options
  and examples.
 }
}
\details{
  The function is a basic implementation of the iterated sure
  independent screening method described
  in Gorst-Rasmussen & Scheike (2011). Briefly, the algorithm does the following:
  \enumerate{
     \item Recruits the \code{nsis} most relevant covariates by ranking them according to the univariate ranking
     method described by \code{rank}.
     \item Selects, using \code{ahazpen} with penalty function described
  in \code{penalty}, a model among the
  top two thirds  of the \code{nsis} most relevant covariates. Call the
  size of this model \eqn{m}.
     \item Recruits `\code{nsis} minus \eqn{m}' new covariates among the non-selected
  covariates by ranking their relevance according to the univariate
  ranking method described in \code{rank}, adjusted for the already
  selected variables (using an unpenalized semiparametric additive
  hazards model).
}
Steps 2-3 are iterated for  \code{maxloop} times, or until \code{nsis} covariates has been recruited, or until the
set of selected covariate is stable between two iterations; whichever
  comes first.

The following choices of ranking method exist:
\itemize{
\item { \code{rank="FAST"} corresponds to ranking, in the initial
  recruitment step only, by the basic FAST- statistic
described in Gorst-Rasmussen & Scheike (2011). If \code{do.isis=TRUE}
then the algorithm sets \code{rank="z"} for subsequent rankings.}
\item{ \code{rank="coef"} corresponds to ranking by absolute value of
(univariate) regression coefficients, obtained via \code{ahaz} 
  }
\item  \code{rank="z"} corresponds to ranking by the \eqn{|Z|}-statistic of
the (univariate) regression coefficients, obtained via \code{ahaz}
\item{ \code{rank="crit"} corresponds to ranking by the size
  of the decrease in
the (univariate) natural loss function used for estimation by \code{ahaz}.
  }
}
}
\value{
  An object with S3 class \code{"ahazisis"}.
  \item{call}{The call that produced this object.}
  \item{initRANKorder}{The initial ranking order.}
  \item{detail.pickind}{List (of length at most \code{maxloop}) listing
    the covariates
    selected in each recruitment step.}
  \item{detail.ISISind}{List (of length at most \code{maxloop}) listing
    the covariates
    selected in each variable selection step.}
  \item{detail.ISIScoef}{List (of length at most \code{maxloop})
    listing the estimated penalized regression coefficients corresponding to
    the indices in \code{detail.ISISind}.}
  \item{SISind}{Indices of covariates selected in the
  initial recruitment step.}
\item{ISISind}{Indices of the final set of covariates selected by the
  iterated algorithm.}
\item{ISIScoef}{Vector of the penalized regression coefficients of the
  covariates in  \code{ISISind}.}
\item{nsis}{The argument \code{nsis}.}
\item{do.isis}{The argument \code{do.isis}.}
\item{maxloop}{The argument \code{maxloop}.}
}
\references{
Gorst-Rasmussen, A. & Scheike, T. H. (2011). \emph{Independent screening for
  single-index hazard rate models with ultra-high dimensional features.}
  Technical report R-2011-06, Department of Mathematical Sciences, Aalborg University.
}
\seealso{\code{\link{print.ahazisis}}, \code{\link{ahazpen}}, \code{\link{ahaz.adjust}}}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Basic ISIS/SIS with a single step
set.seed(10101)
m1 <- ahazisis(surv,X,maxloop=1,rank="coef")
m1
# Indices of the variables from the initial recruitment step
m1$SISind
# Indices of selected variables
m1$ISISind
# Check fit
score <- X[,m1$ISISind]\%*\%m1$ISIScoef
plot(survfit(surv~I(score>median(score))))

}
\keyword{survival}

 
