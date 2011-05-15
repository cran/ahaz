"setupplot"<-function()
  par(cex.main=0.9,tcl=-0.3,mgp=c(1.5,0.5,0),mai=c(0.75,.5,.75,0.1))

"plot.tune.ahazpen" <- function(x, ...)
{
    ## Purpose: Plot 'ahazpencv' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   x     : 'ahazpencv' object
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

  if(x$tune$type=="CV"){
    layout(cbind(1:2), heights=c(8,1))
    setupplot()
    mar<-par()$mar
    plot(1,1,type="n", xlab = expression(lambda),
         ylim = range(c(x$tuneup,x$tunelo)),xlim=rev(range(x$lambda)),
         ylab="Cross-validation score",log = "x", ...)
    error.bars(x$lambda, x$tuneup, x$tunelo, width = 0.001, col = "darkgrey")
    abline(v = x$lambda.min,lty = 2)
    points(x$lambda, x$tunem, pch = 20, col = 2)
    vv<-(c(1,diff(x$df))!=0)
    axis(side = 3, at = x$lambda[vv], labels = paste(x$df[vv]), tick = TRUE, line = 0)
    mtext("# nonzero parameters", side = 3, line = 2, adj = .5, cex = 0)
    par(mar=c(0,mar[2],0,0))
    plot.new()
    legend("topleft",legend=c("Cross-validation score",expression(paste("Optimal ",lambda,sep=""))),
           lty=c(0,2),pch=c(20,-1),col=c("red","black"),bg="white",horiz=TRUE)

  } else {
    layout(cbind(1:2), heights=c(8,1))
    setupplot()
    mar<-par()$mar
    plot(1,1,type="n", xlab = expression(lambda),
         ylim = range(x$tunem),xlim=rev(range(x$lambda)),
         ylab="BIC",log = "x", ...)
    abline(v = x$lambda.min,lty = 2)
    points(x$lambda, x$tunem, pch = 20, col = 2)
    vv<-(c(1,diff(x$df))!=0)
    axis(side = 3, at = x$lambda[vv], labels = paste(x$df[vv]), tick = TRUE, line = 0)
    mtext("# nonzero parameters", side = 3, line = 2, adj = .5, cex = 0)
    par(mar=c(0,mar[2],0,0))
    plot.new()
    legend("topleft",legend=c("BIC value",expression(paste("Optimal ",lambda,sep=""))),
           lty=c(0,2),pch=c(20,-1),col=c("red","black"),bg="white",horiz=TRUE)
           
  }
}
