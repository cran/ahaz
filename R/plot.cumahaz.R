"plot.cumahaz"<-function(x, ...)
  {
    ## Purpose: plot cumulative hazard function
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object: 'cumahaz' object
    ##   ...   : additional arguments to plot function
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
    mint <- min(x$times)
    maxt <- x$times[max(which(x$event != 0))]
    plot(stepfun(x$times, c(0, x$cumhaz)), main = "", xlab = "Time",
         ylab="Cumulative baseline hazard", do.points = FALSE, xlim = c(mint, maxt), ...)
  }
