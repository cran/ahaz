"[.ahazpen"<-function(object,i)
  {
    ## Purpose: Extract ahazpen object at ith lambda-value
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object        : ahazpen object
    ##   i             : lambda index
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    if(missing(i))
      i<-1:length(object$lambda)
    
    small<-object
    small$beta<-object$beta[i,]
    small$lambda<-object$lambda[i]
    small$df<-object$df[i]

    return(small)
  }


"ahazpen" <- function(surv, X, weights,  standardize = TRUE,  penalty=lasso.control(), dfmax = nobs-1,
                      lambda, lambda.minf = ifelse(nobs < nvars,0.01, 1e-4), nlambda = 100,
                      penalty.wgt = NULL, keep = NULL, control=list()){
    ## Purpose: Penalized semiparametric additive hazards regression
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv          : Surv object (right censored/counting process)
    ##   X             : Numeric matrix of covariates
    ##   weights       : Weight vector (nonnegative)
    ##   standardize   : Standardize X?
    ##   penalty       : Control function for initializing penalty
    ##   lambda        : Optional lambda vector
    ##   lambda.minf   : Optional minimum lambda value
    ##   nlambda       : Size of lambda grid
    ##   penalty.wgt   : Optional penalty factor for each covariate
    ##   keep          : Keep covariates in 'keep' (no penalization) 
    ##   dfmax         : Maximal number of covariates to include in path
    ##   control       : List of control parameters for CCD algorithm
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
  
    this.call<-match.call()
    control <- do.call("ahazpen.fit.control",control)


    if(!missing(lambda) && any(lambda <= 0)){stop("'lambda' must be strictly positive")}
    if(nlambda<1){stop("'nlambda' must be strictly positive")}
    if(any(penalty.wgt < 0)){stop("'penalty.wgt' must be nonnegative")}

    # Internal formatting of data
    data <- ahaz.read(surv=surv, X=X, weights=weights, standardize=standardize)
    nobs <-data$nobs
    nvars <- data$nvars

    # Evaluate 'penalty'
    penalty<-eval(penalty)

    if(is.character(penalty)){
      tmp<-c("lasso.control","sscad.control",penalty)[pmatch(tolower(penalty),c("lasso.control","sscad.control"),nomatch=3)]
      penalty<-get(tmp,mode="function")
    }
    if(is.function(penalty)) penalty <- penalty()
    if(is.null(penalty$type)) {
      print(penalty)
      stop("'penalty' not recognized")
    }
    initsol<-penalty$init.sol(surv=surv,X=X,weights=weights)

    
    if (dfmax < 1) { warning("'dfmax' must be strictly positive")}
    if(lambda.minf <= 0){stop("'lambda.minf' must be strictly positive")}
    if(any(keep > nvars | keep <= 0)){stop("'keep' incorrectly specified")}
    if(is.null(penalty.wgt))
      penalty.wgt <- rep(1, nvars)
    penalize<-rep(1,nvars)
    if(!is.null(keep))
      penalize[keep] <- 0
    # Un-penalize variables to be kept - and apply ada.factor
    penalty.wgt<-penalty.wgt*penalize*abs(penalty$ada.wgt(surv,X,weights))
    if(standardize)
      initsol<-initsol*data$standardize$sd
    
    
    if (nvars > nobs && dfmax >= nobs && penalty$alpha == 1.0) {
      warning("'dfmax'<'nobs' required")
      dfmax<-nobs-1
    }
       
     # Get lambda sequence
    uni <-.C("ahd",
             X       = as.double(data$X),
             inout   = as.double(data$inout),
             tdiff   = as.double(data$tdiff),
             iatrisk  = as.double(data$iatrisk),
             deathyn = as.integer(drop(data$death.yn)),
             n       = as.integer(length(data$tdiff)),
             p       = as.integer(data$nvars),
             D       = numeric(data$nvars),
             d       = numeric(data$nvars))
  
    # Penalty function prefactor/initial solution (for sscad penalty)
    prefactor<-ifelse(is.null(penalty$c),ifelse(penalty$type=="sscad",mean(uni$D),1),penalty$c)

    # Get lambda sequence
    lambda <- ahaz.makelambda(lambda, uni$d, lambda.minf, nlambda, penalty.wgt,initsol*prefactor,penalty)
    nlambda <- length(lambda$lambda)

    # Hard-code pmax (maximal number of features that the CCD algorithm ever considers) for simplicity
    pmax = floor(1.2*nobs)

    # Run CCD algorithm
    start<-Sys.time()
    a <- .C("ahpen",
            X          = as.double(t(data$X)),
            inout      = as.double(data$inout),
            wgt        = as.double(data$wgt),
            tdiff      = as.double(data$tdiff),
            times      = as.double(data$tatrisk),
            iatrisk    = as.double(data$iatrisk),
            deathyn    = as.integer(data$death.yn),
            n          = as.integer(length(data$tdiff)),
            p          = as.integer(data$nvars),
            lambda     = as.double(lambda$lambda),
            nlambda    = as.integer(nlambda),
            thresh     = as.double(control$thresh),
            penalty    = as.double(ifelse(penalty.wgt==Inf,-1,penalty.wgt)),
            maxit      = as.integer(control$maxit),
            estims     = as.double(matrix(0, nrow = nlambda, ncol = data$nvars)),
            dfmax      = as.integer(dfmax),
            pmax       = as.integer(pmax),
            d          = as.double(uni$d),
            diagD      = as.double(uni$D),
            lambdaflag = as.integer(nlambda),
            alpha      = as.double(penalty$alpha),
            iter       = as.integer(FALSE),
            initsol    = as.double(initsol),
            a          = as.double(penalty$a), 
            nsteps     = as.integer(penalty$nsteps),
            prefactor  =  as.double(prefactor))
 
    if(a$iter>=control$maxit)
      warning("Estimates may not have converged (try increasing 'maxit'/'thresh' or decreasing 'dfmax')")

    # Return beta estimates on original scale
    if(standardize)
      beta <- t(matrix(a$est,nrow=data$nvars,ncol=nlambda)/data$standardize$sd)
    else
      beta <- t(matrix(a$est,nrow=data$nvars,ncol=nlambda))

    # Number of nonzero covariates
    df <- apply(beta, 1, function(x){sum(x != 0)})
    # Return results for these indices only
    keep <- intersect(1:a$lambdaflag, lambda$keep)
    # If no results returned, 'dfmax'/'pmax'/'nlambda' may conflict 
    if(length(keep)==0)
      stop("no nonzero regression coefficients in fit")
    
    out <- list("call"=this.call, "beta" = drop(beta[keep,]), "lambda" = lambda$lambda[keep],
                "df" = df[keep], "nobs" = data$nobs, "nvars" = data$nvars,
                "surv" = surv,"iter"=a$iter, "penalty.wgt"=penalty.wgt, "penalty"=penalty)
    out$call<-this.call
    
    return(structure(out, class = "ahazpen"))
  }
