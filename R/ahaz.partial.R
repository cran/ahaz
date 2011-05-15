"ahaz.partial" <- function(surv, X, weights, standardize = TRUE, idx)
  {
    ## Purpose: Calculation of columns 'idx' in the matrices D and B used in
    ##          semiparametric additive hazards model
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   surv       : Surv object (right censored/counting process)
    ##   X          : Numeric matrix of covariates
    ##   weights    : Weight vector (nonnegative)
    ##   standardize: Standardize X?
    ##   idx        : Which columns to calculate?
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen

    this.call <- match.call()


    # Internal formatting of data
    tmp <- ahaz.read(surv=surv, X=X, weights=weights, standardize=standardize)

                                        # Check and sort 'idx'; make logical array marking columns in 'idx'
    if(!is.numeric(idx) || length(idx)<1 || length(idx)>tmp$nvars)
      stop("Incorrect 'idx'")
    idx<-sort(as.integer(idx))
    usethis<-rep(-1,tmp$nvars)
    usethis[idx]<-1
    
    # Call C-routine to calculate matrix elements with index < idx 
    a <- .C("ah",
            X       = as.double(tmp$X),
            inout   = as.double(tmp$inout),
            tdiff   = as.double(tmp$tdiff),
            iatrisk = as.double(tmp$iatrisk),
            deathyn = as.integer(drop(tmp$death.yn)),
            n       = as.integer(length(tmp$tdiff)),
            p       = as.integer(tmp$nvars),
            D       = numeric(length(idx)*tmp$nvars),
            d       = numeric(tmp$nvars),
            B       = numeric(length(idx)*tmp$nvars),
            univar  = as.integer(0),
            usethis = as.integer(usethis),
            nuse = as.integer(length(idx)))

    # Call C-routine to calculate maxtrix elements with index > idx
    b <- .C("ah",
            X       = as.double(tmp$X[nrow(tmp$X):1,]),
            inout   = as.double(tmp$inout),
            tdiff   = as.double(tmp$tdiff),
            iatrisk = as.double(tmp$iatrisk),
            deathyn = as.integer(drop(tmp$death.yn)),
            n       = as.integer(length(tmp$tdiff)),
            p       = as.integer(tmp$nvars),
            D       = numeric(length(idx)*tmp$nvars),
            d       = numeric(tmp$nvars),
            B       = numeric(length(idx)*tmp$nvars),
            univar  = as.integer(0),
            usethis = as.integer(rev(usethis)),
            nuse = as.integer(length(idx)))

    # Merge partial D and B matrices
    combmat<-function(x,y,p,idx)
      {
        out<-matrix(0,nrow=length(idx),ncol=p)
        j<-1
        k<-1
        x<-x[x!=0]
        y<-rev(y[y!=0])
        for(i in 1:length(idx)){
          out[i,idx[i]:1]<-x[j:(j+idx[i]-1)]
          out[i,p:idx[i]]<-y[k:(k+p-idx[i])]
          j<-j+idx[i]
          k<-k+(p-idx[i])+1
        }
        return(out)
      }
    D<-as.matrix(combmat(a$D,b$D,tmp$nvars,idx),ncol=length(idx))
    B<-as.matrix(combmat(a$B,b$B,tmp$nvars,idx),ncol=length(idx))
    d<-a$d[idx]

    if (standardize)
      {
        scale <- tmp$standardize$sd[idx] %*% t(tmp$standardize$sd)
        D <- D * scale
        B <- B * scale
        d <- d * tmp$standardize$sd[idx]
      }
    return(list("call"=this.call,"idx"=idx,"nobs"=tmp$nobs,"nvars"=tmp$nvars,"d"=d,"D"=D,"B"=B))
  }
