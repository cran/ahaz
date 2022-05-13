#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define USE_FC_LEN_T
#include <R_ext/BLAS.h>
#include <R.h>

#define EPS (1.0E-6)
		  
// ########
// # BLAS # 
// ########
extern double F77_NAME(ddot)(const int *n, const double *dx, const int *incx,
			    const double *dy, const int *incy);

#ifdef FC_LEN_T
extern void F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
			   const double *alpha, const double *a, const int *lda,
			   const double *x, const int *incx, const double *beta,
                           double *y, const int *incy, size_t);
#else
extern void F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
			   const double *alpha, const double *a, const int *lda,
			   const double *x, const int *incx, const double *beta,
			   double *y, const int *incy);
#endif

double dot(int n,const double * x,const double * y)
{
  // Purpose: Dot product: BLAS
  // ----------------------------------------------------------------------
  int incx = 1;
  int incy = 1;

  return F77_CALL(ddot)(&n, x, &incx, y, &incy);
}
				  
void matvecmul_blas(const double * M, const double *v1, double *v2, int n, int p, char trans) 
{
  // Purpose: Matrix vector multiplication: BLAS
  // ----------------------------------------------------------------------

  double alpha = 1.0;
  double beta = 0.0;
  int incx = 1;
  int incy = 1;
  int nrow = n;
  int ncol = p;

#ifdef FC_LEN_T
  F77_CALL(dgemv)(&trans, &nrow, &ncol, &alpha, M, &nrow, v1, &incx, &beta, v2, &incy, 1);
#else
  F77_CALL(dgemv)(&trans, &nrow, &ncol, &alpha, M, &nrow, v1, &incx, &beta, v2, &incy);
#endif
  }
			  
static inline void matvecmul_nonblas(double ** M, double * v1, double *v2, int n, int p, char trans) 
{
  // Purpose: Matrix vector multiplication: manual
  // ----------------------------------------------------------------------
  // Assigns v2=M*v1 (if trans=='n'), v2=t(M)*v1 otherwise
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  if(trans=='n')
    {
      for (int i = 0; i < n; i++) 
	v2[i] = 0;
      for (int j = 0; j < p; j++) 
	for (int i = 0; i < n; i++) 
	  v2[i] += M[j][i] * v1[j];
       } 
  else {
    for (int j = 0; j < p; j++) 
      v2[j] = dot(n,v1,M[j]);
  }
}    

// ##########################
// # SOFT-THRESHOLDING/SCAD #
// ##########################

static inline double scad(const double prefactor, const double xx, const double l, const double a)
{
 // Purpose: One-step SCAD penalty
  // ----------------------------------------------------------------------
  // Returns lambda * I(x <= l) + I(x > l) * (a*l-x)_+ / (a-1)
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  double x = prefactor*xx;
  if(fabs(x) <= l)
    return(l);
  else {
    double x2 = a * l - fabs(x);
  if(x2 < 0)
    return(0);
  return(l * x2 / ((a - 1) * l) );
  }
}

 static inline double soft_thr(const double x, const double y) {
  // Purpose: Soft thresholding
  // ----------------------------------------------------------------------
  // Returns sign(x)(|x|-y)_+
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen
   return (fabs(x) > y) ? ((x > 0) ? x - y : x + y) : 0;
  }

static inline double signum(const double x) {
  // Purpose: Soft thresholding
  // ----------------------------------------------------------------------
  // Returns sign(x)
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen
  return ( x > 0 ) ? 1 : ( ( x < 0 ) ? -1 : 0);
}


// ###########################
// # SORTING W/PERMUT RETURN #
// ###########################
struct pair {
  double a;
  int b;
};

int compare (const void *const first, const void *const second)
{
  if (((const struct pair *)first)->a < ((const struct pair *)second)->a)
    return 1;
  else if (((const struct pair *)first)->a > ((const struct pair *)second)->a)
    return -1;
  else
    return 0;
}

void sortix(double *x, int *ix, int n)
{
  // Purpose: Sort x w/index return
  // ----------------------------------------------------------------------
  // Sorts x and saves sorting permutation in ix
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  struct pair * ab;
  ab = malloc(n*sizeof(*ab));
  
  for (int i = 0; i < n; i++) {
    ab [i].a = x[i];
    ab [i].b = i;
  }
  qsort (ab, n, sizeof *ab, compare);
  for (int j = 0; j < n; j++) {
      x[j] = ab[j].a;
      ix[j] = ab[j].b;
    }
  free(ab);
} 

// ###########################
// # STANDARDIZATION OF DATA #
// ###########################
void scale (double *X, const double *weights, int n, int p,double *mn, 
	     double *iSd, int standardize)
{
  // Purpose: Center and scale a data matrix
  // ----------------------------------------------------------------------
  // Arguments:
  //   X          : data (design matrix flattened by cols)
  //   weights    : observation weights, length n
  //   n          : nrow(X)
  //   p          : ncol(X)
  //   mn         : placeholder for mean, length n
  //   iSd        : placeholder for inverse sd, length n
  //   standardize: scale? Remember, we always center
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  int npar = p;   // Number of obs.
  int nobs = n;   // Number of obs.
  
  double tmp;
  double invsd;
  double sum = 0;
  double * wgt = Calloc(nobs,double);

  for(int i=0;i<nobs;i++) sum += weights[i];
  for(int i=0;i<nobs;i++) wgt[i] = weights[i]/sum;
  

  if(standardize){
  for (int j=0; j < npar; j++) {
    double msq=0;
    for (int i = 0; i < nobs; i++) {
      tmp = wgt[i] * X[j * nobs + i];
      mn[j] += tmp;
      msq += tmp * X[j * nobs + i];
    }
    tmp = mn[j] * mn[j];
    if (nobs > 1) {
      invsd = 1.0 / sqrt(msq - mn[j] * mn[j]);
    } else {
      invsd = 1;
    }
    iSd[j] = invsd;
    for (int i = 0; i < nobs; i++) {
      X[j * nobs + i] = (X[j * nobs + i] - mn[j]) * invsd;
    }
  }
  } else {
    for (int j=0; j < npar; j++) {
      for (int i = 0; i < nobs; i++) {
	tmp = wgt[i] * X[j * nobs + i];
	mn[j] += tmp;
      }
      for(int i = 0; i < nobs; i++) {
	X[j * nobs + i]=(X[j * nobs + i] - mn[j]);
      }
      iSd[j] = 1.0; 
    }
  }
}

// ######################
// # FORMATTING OF DATA #
// ###################### 
struct survdat {
  double times;
  int ixso;
  int osxi;

  double inout;
  double wgt;
  double tdiff;
  double atrisk;
  double iatrisk;
  double tatrisk;
  int deathyn;
  double totobs;
};

void formatsurvdat(double *time1,double *time2, int *event, double *weights, 
		   struct survdat * s, int n,int rightcens)
{
  // Purpose: Format survival data
  // ----------------------------------------------------------------------
  // Internal formatting of censored/truncated survival times
  // ----------------------------------------------------------------------
  // Arguments:
  //   time1    : entry time (0 if not counting process), length n
  //   time2    : survival time, length n
  //   event    : has subject experienced an event, length n?
  //   weights  : observation weights, length n
  //   s        : placeholder for formatted data
  //   n        : number of observations
  //   rightcens:  Right-censoring (1) or counting process (0) 
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  // Scale weights
  double cumsum = 0;
  double sum = 0;

  // Copy arrays and do sorting
  if (rightcens) {
    for (int i = 0; i < n; i++) sum += weights[i];

    int * ix = Calloc(n,int);
    for(int i = 0; i < n; i++) ix[i] = i;
    
    sortix(time2,ix,n);

    // Construct the relevant variables
    for (int i = 0; i < n; i++) {
	s[i].times = time2[i];
	s[i].ixso = ix[i];
	s[ix[i]].osxi = i;
	s[i].inout = weights[ix[i]];
	s[i].wgt = s[i].inout;
	s[i].tdiff = (i < n-1) ? time2[i] - time2[i+1] : 0.0;
	s[i].deathyn = event[ix[i]];
	cumsum+=s[i].inout; s[i].atrisk = cumsum;
	s[i].iatrisk = s[i].atrisk ? 1.0/s[i].atrisk : 0;
	s[i].tatrisk = time2[i];
      }
    Free(ix);
  } else {
      int nOrig = n / 2;
      for(int i=0;i<nOrig;i++) sum += weights[i];

      int * ix = Calloc(n,int);
      double * tmpTimes = Calloc(n,double);
      int * tmpDeathyn = Calloc(n,int);
      double * tmpInout = Calloc(n,double);
      double * tmpTdiff = Calloc(n,double);

      for (int i = 0; i < nOrig; i++) {
	  tmpTimes[i] = time1[i];
	  tmpDeathyn[i] = 0;
	  tmpInout[i] = -1.0;
	  tmpTdiff[i] = time2[i]-time1[i];
	}
      for (int i = 0; i < nOrig; i++) {
	  tmpTimes[i+nOrig] = time2[i];
	  tmpDeathyn[i+nOrig] = event[i];
	  tmpInout[i+nOrig] = 1.0;
	  tmpTdiff[i+nOrig] = 0.0;
	}

      sortix(tmpTimes,ix,n);

      for(int i = 0; i < n; i++) {
	  s[i].times=tmpTimes[i];
	  s[i].ixso = ix[i];
	  s[ix[i]].osxi = i;

	  s[i].deathyn=tmpDeathyn[ix[i]];
	  s[i].tdiff = (i< n-1) ? tmpTimes[i]-tmpTimes[i+1] : 0;
	  s[i].inout = tmpInout[ix[i]]*weights[ix[i]];
	  s[i].wgt = weights[ix[i]];
	  cumsum += s[i].inout; s[i].atrisk = (i< n-1) ? cumsum : 1.0;
	  s[i].iatrisk = s[i].atrisk ? 1.0 / s[i].atrisk : 0.0;
	  s[i].tatrisk = tmpTdiff[ix[i]]; 
	}

      Free(ix); Free(tmpDeathyn); Free(tmpInout);  Free(tmpTdiff);
      Free(tmpTimes); 
      }
}

// #############################
// # FIT ADDITIVE HAZARD MODEL #
// #############################
void aha (double * X, double *time1, double *time2, int *event, 
	  double *weights, int *n, int *p, double * d, double *D, 
	  double *B, int * getB, int *univariate, int *usethis, int *rightcens)
{
  // Purpose: Fit additive hazards model
  // ----------------------------------------------------------------------
  // Computes D, d and B (optional) for the additive hazards model
  // ----------------------------------------------------------------------
  // Arguments:
  //   X         : Design matrix (flattened by columns)
  //   time1     : Entry time (0 unless counting process data)
  //   time2     : Survival time
  //   event     : Does subject experience an event?
  //   weights   : Observaiton weights
  //   n         : nrow(X)
  //   p         : ncol(X)
  //   d         : Placeholder for d-vector
  //   D         : Placeholder for D-matrix
  //   B         : Placeholder for B-matrix
  //   getB      : Should we calculate B (1) or not (0)
  //   univariate: Return univariate estimes (1) or "full" (0)
  //   usethis   : Vector of indices (C-style) of variables for which to get d/D/B
  //   rightcens : Right-censoring (1) or counting process (0)
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  // If delayed entry, expand data to 2*n obs
  int nobs = *rightcens ? *n : 2*(*n);
  int npar = *p;

  // Format data
  struct survdat * s;
  s = malloc(nobs*sizeof(*s));

  // Do scaling (we must ALWAYS center variables!)
  double * mn = Calloc(npar,double);
  double * iSd = Calloc(npar,double);
  scale(X,weights,nobs,npar,mn,iSd,0);

  // Format survival data
  formatsurvdat(time1,time2,event,weights,s,nobs,*rightcens);

  // Containers for n-vectors for D/B calculations
  double * dotD = Calloc(nobs,double);
  double * dotDunso = Calloc(nobs,double);
  double * dotB = Calloc(nobs,double);
  double * dotBunso = Calloc(nobs,double);
  
  // Container for sorted row of X
  double * XjSo = Calloc(nobs,double);

  // Univariate requested?
  int P = *univariate ? 0 : npar;
  int jj = 0;
  // How many columns to fetch?
  int nuse = 0;
  for (int j = 0; j < npar; j++) {
    if(usethis[j]) {
      nuse++;
    }
  }
  for(int j = 0; j < npar; j++) {
    if (usethis[j]) { 
      double cD = 0; double covsum = 0;  double tmp=0;
      double cB1 = 0; double cB2 = 0;
      d[jj]=0;
      // Sort current row of X
      for (int i = 0; i < nobs; i++) {
	XjSo[i]=X[j * nobs + s[i].ixso];
      }
      for (int i = 0; i < nobs; i++) {
	covsum += XjSo[i] * s[i].inout;
	tmp = covsum * s[i].iatrisk;
	// Get D 
	dotD[i] = s[i].wgt * XjSo[i] * s[i].tatrisk + s[i].inout * cD; 
	cD += tmp * s[i].tdiff;
	D[jj*P+j] += dotD[i] * XjSo[i];
	if( s[i].deathyn )  {
	  d[jj] +=  s[i].inout * (XjSo[i] - tmp);	
	}
	// Get B 
	if (*getB) {
	  dotB[i] =  s[i].inout * (cB1 - cB2);
	  if ( s[i].deathyn )  {
	    dotB[i] += s[i].wgt * (XjSo[i] - tmp);
	    cB1 +=  s[i].wgt * s[i].iatrisk * XjSo[i] ;
	    cB2 += s[i].wgt * tmp * s[i].iatrisk;
	  }
	  B[jj*P+j] += dotB[i] * XjSo[i];
	}
      }
      
      // Get the remaining matrix elements
      if (*univariate == 0) {
	// "Unsort" 
	for (int i = 0; i < nobs; i++) {
	  dotDunso[i] = dotD[s[i].osxi];
	  dotBunso[i] = dotB[s[i].osxi];
	}
	
	// If getting everything, just calculate lower tri
	int jmax = (nuse == npar) ? j : npar;
	for (int k = 0; k < jmax; k++) {
	  D[jj*npar+k] = dot(nobs,dotDunso,X + k * nobs);
	  if(*getB) {
	    B[jj*npar+k] = dot(nobs,dotBunso,X + k * nobs);
	  }
	}
      }
      jj++;
    } 
  }

  Free(dotD); Free(dotDunso); Free(XjSo);
  Free(dotB); Free(dotBunso); Free(iSd); Free(mn);
  free(s);
}




// ################################
// # GET LAMBDA_MAX (SSCAD/LASSO) #
// ################################
double sscadmaxlam(double *initsol, double prefactor, double *d, double a, int p, double *penalty)
{
  // Purpose: Get maximal regularization parameter for sscad/lasso
  // ----------------------------------------------------------------------
  //   initsol  : Initial solution (0 except for SSCAD penalty), length p
  //   prefactor: Prefactor (1 except for SSCAD penalty)
  //   d        : d-vector, length p
  //   a        : Constant a (1 except for SSCAD penalty)
  //   p        : Number of variables
  //   penalty  : Variable weight (differential penalization), length p
  // Author: Anders Gorst-Rasmussen

  double tmp; double sabs; double iabs; 
  double max = 0; double denom;
  for(int j = 0; j < p;j++) {
    if (penalty[j] >= 0) {
      denom = (penalty[j] == 0) ? 1 : penalty[j];
      sabs = fabs(d[j] / denom);
      iabs = fabs(prefactor * initsol[j]);
      tmp = (sabs >= iabs) ? sabs : sabs + (iabs-sabs) / a;
      if (tmp > max) {
	max=tmp;
      }
    }
  }
  return max;
}

// ####################################
// # PENALIZED ADDITIVE HAZARDS MODEL #
// ####################################
void ahapen(double * X, double *time1, double *time2, int *event, 
	    double *weights, int *n,  int *p , int *standardize, 
	    double *lambdaminf, int *nlam, double *lam,  double *thresh, 
	    int *maxit,double *estims, int *dfmax, int *pmax, int *lambdaflag, 
	    int *error, int *activeidx, int *rightcens, double *penalty, double *alpha, 
	    double *initsol, double *a, int *passesleft, int *nsteps, double *prefactor)
{
  // Purpose: Calculate regularization path for penalized additive hazards model
  // ---------------------------------------------------------------------------
  // Arguments:
  //   X          : Design matrix (flattened by columns)
  //   time1      : Entry time (0 unless counting process data), length p
  //   time2      : Survival time, length p
  //   event      : Does subject experience an event, length p?
  //   weights    : Observaiton weights, length p
  //   n          : nrow(X)
  //   p          : ncol(X)
  //   standardize: Should we scale variables when penalizing?
  //   lambdaminf : Minimal lambda as fraction of maximal. If lambdaminf>=1, use 
  //                lam specified by user 
  //   nlam       : Number of lambda-values
  //   lam        : User-specified lambda sequence
  //   thresh     : Threshold for convergence. Stop iterations when relative change in 
  //                penalized loss less than thresh
  //   maxit      : Maximal number of passes over data for *all* lambda values
  //   estims     : Placeholder for coefficients (pmax x nlambda)
  //   dfmax      : Maximal number of variables to include before returning
  //   pmax       : Maximal number of variables to ever consider for inclusion
  //   lambdaflag : For saving number of lambda-values when exiting prematurely
  //   error      : For returning (nonfatal) error messages to R
  //   activeidx  : Index of active variables (combine with estims to get estimates)
  //   rightcens  : Right-censoring (1) or counting process (0)
  //   penalty    : Variable weight (differential penalization), length p
  //   alpha      : alpha for elasticnet penalty
  //   initsol    : Initial solution, length p (0 except for SSCAD penalty)
  //   a          : Constant a (1 except for SSCAD penalty)
  //   passesleft : Maximal number of passes over the data, for all lambda values
  //   nsteps     : Number of steps (1 except for SSCAD)
  //   prefactor  : Prefactor (1 except for SSCAD
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  // If delayed entry, expand data to 2*n obs
  int npar = *p;  
  int nobs = *rightcens ? *n : 2*(*n);


  // ......INITIALIZATION......
 
  // Sort and format data
  struct survdat * s;
  s = malloc(nobs*sizeof(*s));
  formatsurvdat(time1,time2,event,weights,s,nobs,*rightcens);
 
  // Do scaling (we *must* center for our formulas to work!)
  double * mn = Calloc(npar,double);
  double * iSd = Calloc(npar,double);
  scale(X,weights,nobs,npar,mn,iSd,*standardize);

  // Calculate d and find lambda_max
  double * d = Calloc(npar,double);
  for (int j = 0; j < npar; j++) {
    double covsum = 0;
    for (int i = 0; i < nobs; i++) {
      covsum += X[j * nobs + s[i].ixso]  * s[i].inout;
      if (s[i].deathyn) {
	  d[j] += s[i].inout * (X[j * nobs + s[i].ixso] - covsum * s[i].iatrisk);
      }	
    }
  }

  // Maximum lambda - use slightly larger value than necessary for stability
  double maxlam = (1 + EPS) * sscadmaxlam(initsol, *prefactor, d, *a, npar, penalty) / (*alpha);
 

  // If no user-specified lambda sequence, use exponentially decreasing
  int extralam = 0;
  int nlambda;
  double *lambda;
  if (*lambdaminf < 1.0) {
    nlambda = *nlam;
    lambda = Calloc(nlambda,double);
    for (int l = 0; l < nlambda ; l++) {
      lambda[l] = maxlam * pow(*lambdaminf, l / (nlambda - 1.0));
      lam[l] = lambda[l];
    }
  }
  else { // If user-specified lambda, force cold start (10 extra lambda values - better than nothing)
    int *ix = Calloc(*nlam,int);
    sortix(lam,ix,*nlam);
    Free(ix);
    extralam = lam[0] < maxlam ? 10 : 0;
    nlambda = *nlam+extralam;
    lambda = Calloc(nlambda,double);
    
    for (int l = 0; l < nlambda; l++) {
      if(l < extralam){
	lambda[l] = maxlam*pow(lam[0] / maxlam,1.0*l  / (extralam));
      }
      else
	lambda[l] = lam[l-extralam];
    }

  }
    
  // ......CCD PART STARTS HERE......
  double * coefActive = Calloc(npar, double);    // Estimates in active set  
  int * inStrong = Calloc(npar,int);             // Index of feature in strong set (if applicable)
  int * jaTOjc = Calloc(npar,int);               // Active variable index -> complete index
  int * jsTOjc = Calloc(npar,int);               // Strong variable index -> complete index
  int *  inActive =  Calloc(npar,int);           // Logical: is variable in active set?
  int *  isNew =  Calloc(npar,int);              // Logical: is variable new in active set?
  double *  residual = Calloc(nobs,double);      // Residual vector
  double **  D = Calloc(*pmax, double*);         // 'Covariance' matrix
  double * z = Calloc(nobs, double);             // Used for (re)calculating D
  double * riskSc = Calloc(nobs,double);         // Container for risk scores
  double * rs = Calloc(nobs,double);             // For summing p-vectors
  double * res = Calloc(npar,double);            // For summing p-vectors
  double * zz = Calloc(nobs,double);             // For summing p-vectors
  double ** XX = Calloc(*pmax,double*);     // Container for copy of X for active vars

  double * pna = Calloc(npar,double);     // Container for copy of X for active vars
  double * pnb = Calloc(npar,double);     // Container for copy of X for active vars

  int violStrongKKT, violAllKKT;
  double maxdiff, oldLossFun;

  int noActive = 0;
  int sizeStrong = 0;
  double lossFun = 0.0; 
  int noPasses = *maxit;
  for (int j = 0; j < npar; j++) {
    isNew[j] = 1;
  }

  
  // ......LAMBDA LOOP......
  for (int l = 0; l < nlambda; l++) {
    double lambdaomalpha = lambda[l] * (1 - *alpha);
    for (int stp = 0; stp < *nsteps; stp++) {
      for (int j = 0; j < npar; j++) {
	pna[j] = penalty[j] * scad(*prefactor, initsol[j], lambda[l] * (*alpha), *a);
	if (l < nlambda - 1) {
	  pnb[j] = penalty[j] * scad(*prefactor, initsol[j], lambda[l+1] * (*alpha), *a);
	} else {
	  pnb[j] = pna[j];
	}
      }
      // ......FULL SET LOOP......
      violAllKKT = 1;
      while (violAllKKT) {
	violAllKKT = 0;
	
	// ......STRONG SET SET LOOP......
	violStrongKKT = 1;
	while (violStrongKKT) {
	  
	  violStrongKKT = 0;	 
	  
	  // ......CCD LOOP......
	  maxdiff = *thresh+1;
	  while (maxdiff > *thresh && noPasses--) {
	    oldLossFun = lossFun;
	    lossFun = 0.0;
	    for (int ja = 0; ja < noActive; ja++) {
	      int j = jaTOjc[ja];
	      double sum = dot(noActive,coefActive,D[ja]);
	      coefActive[ja] = soft_thr(d[j]  -  sum + coefActive[ja] * D[ja][ja], pna[j]) /  (D[ja][ja] + lambdaomalpha);
	      lossFun += coefActive[ja] * (sum - 2  * d[j] + lambda[l] * signum(coefActive[ja]) + 0.5 * lambdaomalpha * coefActive[ja]);
	    }
	      double denom = oldLossFun ?  oldLossFun : 1.0;
	      maxdiff = fabs((lossFun - oldLossFun) / denom);
	  }
	  if (noPasses < 0){
	      *error = 10;
	      *lambdaflag = l - extralam + 1;
	      goto out;
	  }
	  
	  // ......RESIDUAL......
	    matvecmul_nonblas(XX,coefActive,riskSc,nobs,noActive,'n');

	    double newv = 0; double covsum = 0;
	    for (int i = 0; i < nobs; i++)
	      {
		covsum += riskSc[i] * s[i].inout;
		residual[i] = s[i].wgt * riskSc[i] * s[i].tatrisk + s[i].inout * newv;
		newv += covsum * s[i].iatrisk * s[i].tdiff;
	      }

	    for(int i = 0;i < nobs; i++)
	      rs[i] = residual[s[i].osxi];
      
	    // ......KKT CHECK/STRONG SET......
	    for (int jc = 0; jc < npar; jc++) {
	      if (inStrong[jc]) {
		  double tmp = 0;
		  tmp = dot(nobs,rs,X+jc*nobs);

		  // KKT check
		  inActive[jc] = fabs(d[jc] - tmp) >= pna[jc];
		  
		  // If variable becomes active...
		  if (inActive[jc] && isNew[jc]) {
		      // Break out if needed
		      if( noActive >= *pmax){
			*error = 20;
			*lambdaflag = l-extralam+1;
			goto out;  
		      }   

		      // Add new row
		      D[noActive] = Calloc(*pmax, double); 
		      double newv = 0; double covsum = 0;
		      for (int i = 0; i < nobs; i++)
			{
			  covsum += X[jc * nobs + s[i].ixso] * s[i].inout;
			  z[i] = s[i].wgt * X[jc * nobs + s[i].ixso] * s[i].tatrisk + s[i].inout * newv;
			  newv += covsum * s[i].iatrisk * s[i].tdiff;
			}

		      // Copy 'active covariates'
		      XX[noActive]=Calloc(nobs,double);
		      for (int i = 0; i < nobs; i++) {
			XX[noActive][i]=X[nobs * jc+s[i].ixso];
		      }

		      matvecmul_nonblas(XX,z,D[noActive],nobs,noActive+1,'t');

		      // Update remaining rows by transposition
		      for(int ja = 0; ja < noActive; ja++) {
			D[ja][noActive] = D[noActive][ja];
		      }
		      
		      // Update other stuff
		      jaTOjc[noActive] = jc;
		      isNew[jc] = 0;
		      noActive++;
		      violStrongKKT++;
		    }
		}
	    }
	}
	
	// ......KKT CHECK/FULL SET......
	matvecmul_blas(X,rs,res,nobs,npar,'t');
	for (int jc = 0; jc < npar; jc++) {
	  // If variable NOT in strong set - check KKT/strong
	  if(inStrong[jc] == 0 && penalty[jc] >= 0)
	    {
	      double chk = fabs(d[jc] - res[jc]);
	      if (chk  >= 2 * pnb[jc] - pna[jc]) {
		inStrong[jc]++;
		sizeStrong++;
	      if(chk >= pna[jc]){
		violAllKKT++;
	      }
	      }
	    }
	}
      }
      // If multistep, update initsol
      if (*nsteps > 1) {
	for (int jc = 0; jc < npar; jc++) {
	  initsol[jc]=0;
	}	
	for(int ja = 0; ja < noActive; ja++)
	  initsol[jaTOjc[ja]] = coefActive[ja];
      }
    }
      
    // ......SAVE RESULTS ......
    int noNonzero = 0;
    for (int ja = 0; ja < noActive; ja++) {
      noNonzero += (coefActive[ja]!=0);
      if(l>=extralam) {
	// Un-standardize if needed
	estims[(l-extralam) * (*pmax) + ja] = coefActive[ja] * iSd[jaTOjc[ja]];
      }
    }	
    
    if (noNonzero > *dfmax && *lambdaminf < 1) {
      *lambdaflag = l-extralam+1;
      goto out;
    } else if(noNonzero > nobs - 1 && *alpha >= 1.0 && *lambdaminf < 1) {
      *error = 30;
      *lambdaflag = l - extralam + 1;
      goto out;
    }
  }
    
 out:;
  
  *passesleft = noPasses+1;
  // Return indices of active vars
  for (int ja = 0; ja < noActive; ja++) {
    activeidx[ja] = jaTOjc[ja] + 1;
  }
  
  
  for (int j = 0; j < noActive; j++) {
    Free(D[j]);
    Free(XX[j]);
  }
  Free(D);
  Free(XX);
  
  Free(coefActive); Free(inStrong); Free(jaTOjc); Free(jsTOjc);
  Free(inActive); Free(isNew); Free(residual); Free(z); Free(riskSc);
  Free(rs); Free(res); Free(zz); Free(d); 
  Free(lambda); Free(mn); Free(iSd); Free(pna);Free(pnb);
  
  free(s);    
}


// #############
// # RESIDUALS #
// #############
void ahresid (double *start, double *end, double *status, double *X, 
	      double *Zbar, double *times, double *tdiff, double *breslow, 
	      double *beta, int *ntimes,int *p, int *nobs, double *resid, double *wgt)
{
  // Purpose: Get integrated martingale residuals
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen 

  double tmp;
  for(int t = 0; t < *ntimes - 1; t++)
    {
      for(int i=0; i < *nobs; i++)
	{
	  if(end[i] >= times[t] && start[i] <= times[t+1])
	    {
	      if(status[i] == 1 && times[t] == end[i])
		{
		  for(int k = 0; k < *p; k++)
		    resid[i*(*p)+k] +=  wgt[i] *  (X[i*(*p)+k] - Zbar[t*(*p)+k]);
		}
	      tmp = 0;
	      for(int k = 0; k < *p; k++)
		tmp += X[i * (*p) + k] * beta[k];
	      for(int k = 0; k < *p; k++)
		resid[i * (*p)+k] +=  wgt[i] * (Zbar[t * (*p)  +k] - X[i* (*p) + k]) * (breslow[t] + tmp * tdiff[t]);
	    }
	}
    }
}

// #####################
// # BRESLOW ESTIMATOR #
// #####################
void ahbreslow  (double *X, double *tdiff,double *inout,double *iatrisk, int *deathyn, 
	       int *n,int *p, double *beta,double *bresl,double *zbar)
{ 
  // Purpose: Breslow estimate of cumulative hazard
  // ----------------------------------------------------------------------
  // Arguments:
  //   X      : data (design matrix flattened by rows)
  //   inout  : indicator, entry (+obs. weight)/exit (-obs. weight)
  //   tdiff  : difference between event times
  //   iatrisk: 1/atrisk
  //   deathyn: is event time a death time?
  //   n      : number of event times
  //   p      : number of vars.
  //   S      : placeholder for S
  //   s      : --"-- s
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  double *covsum;
  covsum= Calloc(*p,double);
 
  for(int i = 0; i < *n; i++)
    {
      if(deathyn[i])
	bresl[i] += inout[i] * iatrisk[i];
 
      for(int k = 0; k < *p; k++)
	{
	  covsum[k] += inout[i] * X[i * (*p) + k];
	  zbar[i * (*p) + k] = covsum[k]* iatrisk[i];
	  bresl[i] -= zbar[i * (*p)+k] * tdiff[i] * beta[k];
	}     
    }  
  Free(covsum);
}
