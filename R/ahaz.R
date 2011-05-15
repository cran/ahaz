"ahaz.buildmat" <- function(x, p)
  {
    ## Purpose: Symmetric matrix from 1D flipped upper triangle
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x: Numeric p*(p+1)/2 vector (backwards rowwise matrix entries)
    ##   p: Dimension of matrix
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    f <- function(y) { c(rev(y[y != 0]),y[y == 0]) }
    M <- matrix(0, nrow = p, ncol = p)
    M[(!lower.tri(M))] <- x
    M <- apply(M, 2, f)
    out<-M + t(upper.tri(M) * M)

    return(out)
  }

"ahaz" <- function(surv, X, weights, standardize = TRUE, univariate = FALSE, robust=FALSE)
  {
    ## Purpose: Semiparametric additive hazards regression (Lin & Ying 1994)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   standardize: Standardize X?
    ##   univariate : Do marginal analyses?
    ##   robust     : Setup robust estimation of covariances
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    this.call <- match.call()

    if(univariate && robust)
      stop("option 'robust' not supported for univariate regressions")

    # Validity checks/internal formatting
    tmp <- ahaz.read(surv, X, weights, standardize)

    # Fit additive hazard model
    a <- .C("ah",
            X       = as.double(tmp$X),
            inout   = as.double(tmp$inout),
            tdiff   = as.double(tmp$tdiff),
            iatrisk  = as.double(tmp$iatrisk),
            deathyn = as.integer(drop(tmp$death.yn)),
            n       = as.integer(length(tmp$tdiff)),
            p       = as.integer(tmp$nvars),
            D       = numeric(ifelse(univariate,tmp$nvars,tmp$nvars*(tmp$nvars+1)/2)),
            d       = numeric(tmp$nvars),
            B       = numeric(ifelse(univariate,tmp$nvars,tmp$nvars*(tmp$nvars+1)/2)),
            univar  = as.integer(univariate),
            usethis = as.integer(rep(1,tmp$nvars)),
            nuse = as.integer(ceiling((tmp$nvars+1)/2)))

    # Re-scale estimates if 'standardize=TRUE'
    if (standardize)
      {
        tmp$X <- sweep( tmp$X, 1, tmp$standardize$sd, '*')
        tmp$X <- sweep( tmp$X, 1, tmp$standardize$mean, '+')
        if (!univariate) {
          scale <- tmp$standardize$sd %*% t(tmp$standardize$sd)
          D <- ahaz.buildmat(a$D, tmp$nvars) * scale
          B <- ahaz.buildmat(a$B, tmp$nvars) * scale
          d <- a$d * tmp$standardize$sd
        } else {
          D <- a$D * tmp$standardize$sd^2
          B <- a$B * tmp$standardize$sd^2
          d <- a$d * tmp$standardize$sd
        }
      } else {
            if (!univariate) {
              D <- ahaz.buildmat(a$D, tmp$nvars)
              B <- ahaz.buildmat(a$B, tmp$nvars)
              d <- a$d
            } else {
              D <- a$D
              B <- a$B
              d <- a$d
            }
          }
    
    names(d)<-tmp$colnames
    if(!univariate)
      {
        colnames(D)<-rownames(D)<-tmp$colnames
        colnames(B)<-rownames(B)<-tmp$colnames
      } else {
        names(D)<-names(B)<-tmp$colnames
      }
    out <-structure(list("call" = this.call, "nobs" = tmp$nobs, "nvars" = tmp$nvars,
                "D" = D, "d" = d, "B" = B, "univar" = univariate,
               "data" = tmp, "robust" = 0), class="ahaz")

    # Calculate robust variance estimator?
    if (robust)
      {
        M <- predict(out, type = "residuals")
        B <- t(M) %*% M
        out$B <- drop(B)
        out$robust <- TRUE
      }
    
    return(out)
  }
