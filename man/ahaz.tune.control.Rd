\name{ahaz.tune.control}
\alias{ahaz.tune.control}
\alias{cv.control}
\alias{bic.control}
\title{Tuning controls for regularization}
\description{Define the type of tuning method used for regularization. Currently only used by \code{tune.ahazpen}. }
\usage{
# Cross-validation
cv.control(nfolds=5, reps=1, foldid=NULL, trace=FALSE)

# BIC-inspired
bic.control(factor = function(nobs){log(nobs)})
}
\arguments{
  \item{nfolds}{Number of folds for cross-validation. Default is
    \code{nfolds=5}. Each fold must
  have size > 1, i.e. \code{nfolds} must be
  less than half the sample size.}
  \item{reps}{Number of repetitions of cross-validation with
    \code{nfolds} folds. Default is \code{rep=1}. A \code{rep}
    larger than 1 can be useful to reduce variance of cross-validation
  scores.}
  \item{foldid}{An optional vector of values between 1 and \code{nfolds}
    identifying the fold to which each observation belongs. Supercedes
    \code{nfolds} and \code{rep} if supplied.}
  \item{trace}{Print progress of cross-validation. Default is \code{trace=FALSE}.}
 \item{factor}{Defines how strongly the number of nonzero penalty
   parameters penalizes the score in a BIC-type criterion; see the details.}
}

\details{
  For examples of usage, see \code{\link{tune.ahazpen}}.
  
  The regression coefficients of the semiparametric additive hazards
  model are estimated by solving a linear system of estimating equations of the form
  \eqn{D\beta=d}{D*beta = d} with respect to \eqn{\beta}{beta}. The natural loss function
  for such a linear function is of the least-squares type
  \deqn{L(\beta)=\beta' D \beta -2d'\beta.}{L(beta) = beta' * D * beta -
  2 * d' * beta}
  This loss function is used for cross-validation as described by Martinussen
  & Scheike (2008).

  Penalty parameter selection via a BIC-inspired approach was described by
   Gorst-Rasmussen & Scheike (2011). With \eqn{df} is the degrees of freedom and \eqn{n} the number of
  observations, we consider a BIC inspired criterion of the form
  \deqn{BIC = \kappa L(\beta) + df\cdot factor(n)}{BIC = kappa * L(beta)
   + df * factor(n) }
 where \eqn{\kappa}{kappa} is a scaling constant included to remove dependency on the
   time scale and better mimick the behavior of a `real' (likelihood) BIC. The default \code{factor=function(n){log(n)}} has
    desirable theoretical properties but may be conservative in practice.
}


\value{An object with S3 class \code{"ahaz.tune.control"}. 
 \item{type}{Type of penalty.}
 \item{factor}{Function specified by \code{factor}, if applicable}
 \item{getfolds}{A function specifying how folds are calculated, if applicable.}
 \item{rep}{How many repetitions of cross-validation, if applicable.}
 \item{trace}{Print out progress?}
}

\references{ Gorst-Rasmussen, A. & Scheike,
T. H. (2011). \emph{Independent screening for single-index hazard rate
models with ultra-high dimensional features.}  Technical report
R-2011-06, Department of Mathematical Sciences, Aalborg University.

  Martinussen, T. & Scheike, T. H. (2008). \emph{Covariate selection for
  the semiparametric additive risk model.} Scandinavian Journal of
Statistics; \bold{36}:602-619.  }
\seealso{\code{\link{tune.ahazpen}}}

\keyword{misc}

 
