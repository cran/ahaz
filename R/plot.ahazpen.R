"plot.ahazpen"<-function(x, scale = c("coef","lambda"), labels = FALSE, df = TRUE, ...)
  {
    ## Purpose: plot regularization path from 'ahazpen' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x    : 'ahazpen' object
    ##   type : 'lambda' - first axis is lambda (log scale)
    ##           'coefficients' - first axis is L1 norm
    ##   label: show beta indices in right margin
    ##   df   : show degrees of freedom in top margin 
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    scale <- match.arg(scale)
    if (length(x$lambda) == 1)
      stop("No coefficient paths to plot!")

    beta<-as.matrix(x$beta)
   
    include <- apply(beta, 2, function(x){any(abs(x) > .Machine$double.eps)})
    toplot <- apply(beta, 1, function(x){any(abs(x) > .Machine$double.eps)})

    if(!sum(include))
      stop("All coefficient paths identically zero!")
    
    col <- 1:sum(include)
    if (scale == "lambda"){
      z <- x$lambda
      matplot(z, beta[,include],type = "l", xlab = expression(lambda),
            ylab = "Regression coefficients",lty = 1,col = col,pch = 1,xlim=rev(range(z)),log = "x",...)
    } else {
      z <- apply(as.matrix(beta[,include]), 1, function(x){sum(abs(x))})
      matplot(z,beta[,include], type = "l", xlab = "L1 norm",
              ylab = "Regression coefficients", lty = 1, col = col, pch = 1,...)
    }
    if(labels)
      axis(4, at = beta[max(which(toplot)), include], labels = as.character(1:sum(include)))
    if(df)
      {
        mtext("Number of nonzero coefficients", side = 3, line = 2, adj = .5, cex = 0)
        vv<-(c(1,diff(x$df))!=0)
        axis(side = 3, at = z[vv], labels = paste(x$df[vv]), tick = TRUE, line = 0)
      }
  }
