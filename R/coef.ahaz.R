"coef.ahaz"<-function(object, ...)
  {
    ## Purpose: beta estimates from 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object  : 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    return(predict(object, type="coef"))
}
