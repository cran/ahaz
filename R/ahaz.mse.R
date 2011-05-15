"ahaz.mse"<-function(x, beta)
  {
    ## Purpose: MSE from ahaz object and given beta
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahaz' object
    ##   beta  : beta estimate
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    return((t(beta) %*% x$D %*% beta - 2 * t(beta) %*% x$d))
  }
