"ahazpartial" <- function(surv, X, weights, standardize = TRUE, idx)
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
            S       = numeric(length(idx)*tmp$nvars),
            s       = numeric(tmp$nvars),
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
            S       = numeric(length(idx)*tmp$nvars),
            s       = numeric(tmp$nvars),
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
    S<-combmat(a$S,b$S,tmp$nvars,idx)
    B<-combmat(a$B,b$B,tmp$nvars,idx)
    
    return(list("S"=S,"B"=B,"call"=this.call))
  }



"adjustahaz"<-function(surv,X,idx,method=c("fast","coeff","t","criterion"),standardize=TRUE)
  {
    if(missing(idx) || any(idx<1) || any(idx>ncol(X)))
      stop("'idx' incorrectly specified")
    method<-match.arg(method)
    ipr<-function(x,B=NULL)
      {
        if(is.null(B))
          return(apply(x*x,2,sum))
        else
          return(apply(x*(B%*%x),2,sum))
      }
    
    idx<-sort(idx)
    if(method=="t" || method=="fast")
      { 
        m<-ahazpartial(surv,X,standardize=standardize,idx=idx)
        Sbig<-m$S
        Bbig<-m$B

        morig<-ahaz(surv,X[,idx],standardize=standardize)
        Borig<-morig$B
        Sorig<-morig$S
        
        muni<-ahaz(surv,X,standardize=standardize,univariate=TRUE)
        s<-muni$s
        S<-muni$S

        # Nominator
        iA<-solve(Sorig)
        n<-ncol(iA)
        vf<-iA%*%Sbig
        ik<-1/(S-apply(Sbig*vf,2,sum))
        x1<--ik*t(vf)
        nomin<-x1%*%s[idx]+ik*s

        # Denominator
        denomin<-(ik^2*muni$B+2*ik*apply(Bbig*t(x1),2,sum)+ipr(t(x1),Borig))/nrow(X)
        denomin[idx]<-NA
        nomin[idx]<-NA
        return(abs(nomin/sqrt(denomin)))
        
      } else  {
        m1<-ahazpartial(surv,X,standardize=standardize,idx=idx)
        m2<-ahaz(surv,X[,idx],standardize=standardize)
        m3<-ahaz(surv,X,standardize=standardize,univariate=TRUE)
        SS<-m1$S
        iA<-solve(m2$S)
        S<-m3$S
        s<-m3$s
     
        n<-ncol(iA)
        vf<-iA%*%SS
        ik<-1/(S-apply(SS*vf,2,sum))
        x1<--ik*t(vf)%*%s[idx]
        if (method=="coeff"){
          x2<-ik*s
          na<-x1+x2
          na[idx]<-NA
          return(na)
        }
        else {
          na<--(ik*s^2+2*s*x1+x1^2/ik)
          na[idx]<-NA
          return(na)
        } 
      } 
  }
