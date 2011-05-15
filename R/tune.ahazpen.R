"tune.ahazpen"<-function(surv, X, weights, standardize = TRUE, penalty=lasso.control(), tune=cv.control(),dfmax=nobs-1,lambda,...)
  {
    ## Purpose: Tuning parameter selection for additive hazards lasso
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   standardize: Standardize X?
    ##   K          : Number of folds
    ##   trace      : Print out progress?
    ##   allfolds   : Optional user-specified list of folds
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    this.call <- match.call()
    if(!missing(lambda))
      if(length(lambda)<2)
        stop("'lambda' should have length > 1")

    tune<-eval(tune)
    # Get tuning controls
    if(is.character(tune)){
      tmp<-c("cv.control","bic.control",tune)[pmatch(tolower(tune),c("cv.control","bic.control"),nomatch=3)]
      tune<-get(tmp,mode="function")
    }
    if(is.function(tune)) tune <- tune()
    if(is.null(tune$type)) {
      print(tune)
      stop("'tune' not recognized")
    }
    
    
    # Preliminary 'ahazpen' fit
    nobs<-nrow(X)
    if(missing(lambda))
      fit <-ahazpen(surv=surv,X=X,standardize=standardize,penalty=penalty,dfmax=dfmax,...)
    else
       fit <-ahazpen(surv=surv,X=X,standardize=standardize,penalty=penalty,dfmax=dfmax,lambda=lambda,...)
    nobs<-fit$nobs
    df <- fit$df
    # If no nonzero coefficients in fit, no point in continuing
    if(max(df)==0)
      stop("no nonzero coefficients")

    lambdanew <- fit$lambda
    nlambda <- length(lambdanew)

    # CROSS-VALIDATION
    if(tune$type=="CV")
      {
        foldsused<-list()
        k<-1
        for(j in 1:tune$rep)
          {
            if(tune$trace)
              cat(paste("Repetition: ",j,"/",tune$rep,"\n",sep=""))
            cvf<-tune$getfolds(nrow(X))
            foldsused[[k]]<-cvf;k<-k+1
            if(j==1){
              error <- matrix(0, nrow = cvf$nfolds, ncol = length(lambdanew))
              ll<-matrix(0,nrow=cvf$nfolds,ncol=tune$rep)
            }
            # Slightly increase dfmax to improve chances that we achieve dfmax in CV
            dfmaxnew<-min(floor(max(df)*1.2),nrow(X)-length(cvf$all.folds[[1]])-1)
            for(i in 1:cvf$nfolds)
              {
                if(tune$trace)
                  cat(paste("  Fold: ",i,"/",cvf$nfolds,"\n",sep=""))
                omit <- cvf$all.folds[[i]]
                tmp <- ahazpen(surv=surv[-omit,], X=X[-omit,],lambda=lambdanew,penalty=penalty,dfmax=dfmaxnew,...)  
                ll[i,j]<-length(tmp$lambda)
                beta<-rbind(coef(tmp))
                nzero <- apply(beta, 2, function(x){any(x != 0)})
                test <- ahaz(surv[omit,], X[omit,nzero])
                error[i,1:ll[i,j]] <-error[i,1:ll[i,j]]+apply(beta,1,function(x){ahaz.mse(test, x[nzero])})
              }
          }
        ll<-min(ll)

        tunem  <- apply(cbind(error[,1:ll])/tune$rep,2,mean)
        tunesd <- apply(cbind(error[,1:ll])/tune$rep,2,sd)
        lambda.min <- max(lambdanew[tunem<=min(tunem)])
        
        out <- list("lambda" = lambdanew[1:ll], "tunem" = tunem,
                    "tunesd" = tunesd, "tuneup" = tunem + tunesd,
                    "tunelo" = tunem - tunesd, "lambda.min" = lambda.min,
                    "df" = df[1:ll], "call" = this.call,
                    "tune"=tune,"penalty"=penalty,"foldsused"=foldsused,"nfolds"=cvf$nfolds)
        
        class(out) <- "tune.ahazpen"
        return(out)
      } else if(tune$type=="BIC")
        {
          nz<-apply(coef(fit)!=0,2,any)
          m<-ahaz(surv,X[,nz])
          factor<-tune$factor(nobs)
          invB<-ahaz.ginv(m$B)
          invD<-ahaz.ginv(m$D)
          kappa<-drop(t(m$d)%*%invB%*%m$d/(t(m$d)%*%invD%*%(m$d)))
          bic<-kappa*apply(fit$beta,1,function(x){ahaz.mse(m,x[nz])})+fit$df*factor/nobs
          lambda.min <- max(lambdanew[bic<=min(bic)])
          out<-list("df"=fit$df,"lambda"=lambdanew,"tunem"=bic, "lambda.min"=lambdanew[which.min(bic)],
                     "call" = this.call,"tune"=tune,"penalty"=penalty)
          class(out) <- "tune.ahazpen"
          return(out)
        }  else 
        {
          stop("invalid 'tune'")
        }
  }
