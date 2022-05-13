\name{ahazpen.pen.control}
\alias{ahazpen.pen.control}
\alias{lasso.control}
\alias{sscad.control}
\title{Penalty controls for ahazpen}
\description{Describe the penalty function to be used in the penalized
  semiparametric additive hazards model. Typically only used in a call
  to \code{ahazpen} or \code{tune.ahazpen}.}
\usage{
# (Adaptive) lasso/elasticnet
lasso.control(alpha=1, ada.wgt=NULL)

# Stepwise SCAD
sscad.control(a=3.7, nsteps=1, init.sol=NULL, c=NULL)
}
\arguments{
  \item{alpha}{Elasticnet penalty parameter with default
  \code{alpha=1} corresponding to the standard lasso; see details.}
  \item{ada.wgt}{Optional covariate weights used for fitting the adaptive
  lasso. Default is \emph{not} to use weights, i.e. fit the standard lasso. A
  user-specified \code{init.sol}
  can be a nonnegative vector of length corresponding to the number of covariates in the
  model. \cr For advanced use it may also be specified as a function with
  arguments  \code{surv}, \code{X} and \code{weights} precisely;
  see the details.}
  \item{a}{Parameter of the stepwise SCAD penalty, see details. Default
    is \code{a=3.7}}
 \item{nsteps}{Number of steps in stepwise SCAD. Default is
  \code{nsteps=1}.}
   \item{init.sol}{Optional initial solution for stepwise SCAD
     consisting of  a numerical vector of length corresponding to the number of covariates in the
  model. Default
     is a vector of regression coefficients obtained from
  \code{ahaz} if there are more observations than covariates, zero
  otherwise. 
  For advanced use, \code{initsol} it can also be specified as a \bold{function} with
  arguments  \code{surv}, \code{X} and \code{weights} precisely;
  see the details.}
\item{c}{Optional scaling factor for stepwise SCAD. Usually it is not
  necessary to change supply this; see the details.}
}
\details{
  The lasso/elasticnet penalty function takes the form
  \deqn{p_\lambda(\beta)=\lambda((1-\alpha)\|\beta\|_2 +
    \alpha\|\beta\|_1)}{p(lambda, beta) = lambda*[(1 - alpha)*|beta|_2 + alpha*|beta|_1]}
  where \eqn{0 <\alpha \leq 1}{0 < alpha <= 1}. Choosing \eqn{\alpha<1}{alpha<1} encourages joint selection of correlated covariates and may
  improve the fit when there are substantial correlations among
  covariates. 

  The stepwise SCAD penalty function takes the form \deqn{
  p_\lambda(\beta)=w_\lambda(c|b_1|)|\beta_1|+\ldots+
  w_\lambda(c|b_{nvars}|)|\beta_{nvars}|}{p(lambda, beta) =
  w(lambda, c*|b_1|)
  * |beta_1| + ... + w(lambda, c*|b_nobs|)
  * |beta_nobs|}
where \eqn{b} is some
  initial estimate, \eqn{c} is a scaling factor, and for
  \eqn{I} the indicator function
\deqn{w_\lambda(x)=\lambda I(x \le \lambda) + \frac{(a\lambda
   - x)_+}{a - 1}I(x>\lambda)}{w(lambda, x) = lambda * I(x <= lambda) + I(x > lambda) *
   (a*lambda-x)_+ / (a-1)}
 The scaling factor \eqn{c}  controls how `different' the
 stepwise SCAD penalty is from the standard lasso penalty (and is
  also used to remove dependency of the penalty on the scaling of the time
  axis).

  The one-step SCAD method of Zou & Li (2008) corresponds to taking
  \eqn{b} equal to the estimator derived from \code{\link{ahaz}}. See
  Gorst-Rasmussen & Scheike (2011) for details. By iterating such
  one-step SCAD and updating the initial solution \eqn{b} accordingly,
  the algorithm approximates the solution obtained using full SCAD. Note
  that calculating the full SCAD solution generally leads to a nonconvex
  optimization problem: multiple solutions and erratic behavior of
  solution paths can be an issue.

  The arguments \code{ada.wgt} and \code{init.sol} can be specified as
  functions of the observations. This is convenient, for example, when
  using cross-validation for tuning parameter selection. Such a function
  \strong{must} be specified precisely with the arguments \code{surv},
  \code{X} and \code{weights} and \strong{must} output a numeric vector
  of length corresponding to the number of covariates. \code{ahazpen}
  will take care of scaling so the function should produce output on the
  original scale. See the examples here as well as the examples for
  \code{\link{tune.ahazpen}} for usage of this feature in practice.
  
}



\value{An object with S3 class \code{"ahaz.pen.control"}. 
 \item{type}{Type of penalty.}
 \item{init.sol}{Function specifying the initial solution, if applicable.}
 \item{alpha}{Value of \code{alpha}, if applicable.}
 \item{nsteps}{Number of steps for stepwise SCAD penalty, if
   applicable.}
 \item{a}{Parameter for stepwise SCAD penalty, if applicable.}
  \item{c}{Scaling factor for stepwise SCAD penalty, if applicable.}
 \item{ada.wgt}{Function specifying the weights for the adaptive lasso penalty, if applicable.}
}
\references{
   Gorst-Rasmussen, A. & Scheike, T. H. (2011). \emph{Independent screening for
  single-index hazard rate models with ultra-high dimensional features.}
  Technical report R-2011-06, Department of Mathematical Sciences, Aalborg
  University.
  
    Leng, C. & Ma, S. (2007). \emph{Path consistent model selection in
 additive risk model via Lasso}. Statistics in Medicine; \bold{26}:3753-3770.

  Martinussen, T. & Scheike, T. H. (2008). \emph{Covariate selection
  for the semiparametric additive risk model.} Scandinavian Journal of
Statistics; \bold{36}:602-619.
  
Zou, H. & Li, R. (2008). \emph{One-step sparse estimates in nonconcave
 penalized likelihood models}, Annals of Statistics; \bold{36}:1509-1533.
  }

\seealso{\code{\link{ahazpen}}, \code{\link{tune.ahazpen}}}
\examples{
data(sorlie)

# Break ties
set.seed(10101)
time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Survival data + covariates
surv <- Surv(time,sorlie$status)
X <- as.matrix(sorlie[,3:ncol(sorlie)])

# Fit additive hazards regression model with elasticnet penalty
model <- ahazpen(surv,X,penalty=lasso.control(alpha=0.1),dfmax=30)
plot(model)

# Adaptive lasso with weights 1/|beta_i|^0.5. Note that, although
# we do not use 'weights', it MUST be included as an argument
adafun <- function(surv,X,weights)
 return(1/abs(coef(ahaz(surv,X)))^.5)
model <- ahazpen(surv,X[,1:50],penalty=lasso.control(ada.wgt=adafun))
plot(model)

# One-step SCAD with initial solution derived from univariate regressions
scadfun <- function(surv,X,weights){
 fit <- ahaz(surv,X,univariate=TRUE)
 return(coef(fit))
}
set.seed(10101)
model.ssc <- tune.ahazpen(surv,X,dfmax=30,penalty=sscad.control(init.sol=scadfun))
plot(model.ssc)
}
\keyword{misc}