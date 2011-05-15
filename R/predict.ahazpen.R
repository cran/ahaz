"predict.ahazpen" <- function(object, newX, type = c("coef","lp","residuals","cumhaz"), index = NULL, ...)
  {
    ## Purpose: Prediction for 'ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahazpen' object
    ##   newX  : Covariate values at which to make predictions - required!
    ##   type  : 'linear' risk score
    ##           'residuals' martingale residuals from pen. beta estimate
    ##           'cumhaz' Breslow estimate of cum. haz. from pen. beta estimate
    ##           'mse' MSE from pen. beta estimate
    ##   index : For which lambda index should predictions be made?
    ##           Required for all other types than 'coefficients','lp','mse'
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    type <- match.arg(type)
    if(!is.null(index) && (any(index<1) || !is.null(index) && any(index>length(object$lambda))))
      stop("argument 'index' out of bounds")

    if (type == "coef") {
      if(is.null(index))
        return(object$beta)
      return(object$beta[index,])
    } 
    
    if (missing(newX)) stop("no new data provided")
    if (!is.numeric(newX)) stop("argument 'newX' must be a numeric matrix")
    if(dim(newX)[2]!=object$nvars)
      stop("incorrect dimensions of argument 'newX'")
    if(!is.null(index) && (any(index<1) || !is.null(index) && any(index>length(object$lambda))))
      stop("argument 'index' out of bounds")

    newX<-as.matrix(newX)
    
    if (type=="lp") {
      if(is.null(index))
        return(drop(drop(object$beta)%*%t(newX)))
      return(drop(drop(object$beta)%*%t(newX))[,index])

    }  else {
        if(is.null(index) && length(object$lambda)>1)
          stop("missing argument 'index'")
        if(length(index)>1)
          stop("argument 'index' must have length 1 for option 'type=\"residuals\"'")
        if(dim(newX)[1]!=object$nobs)
          stop("incorrect dimensions of argument 'newX'")
        if(is.null(index) && length(object$lambda)==1)
          index <- 1

        if(is.matrix(object$beta))
          beta <- object$beta[index,]
        else
          beta <- object$beta
        include <- (beta!=0)
        
        if(!is.null(colnames(newX)))
          names <- colnames(newX)[include]
        else
          names <- which(include)

        if(any(include))
          return(ahaz.predbackend(ahaz.read(object$surv, newX[,include], standardize = FALSE),
                                  beta = beta[include], type = type, colnames = names))
        stop("no variables included!")
      }
  }
