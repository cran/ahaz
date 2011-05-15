#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

static inline double scad(const double x, const double l, const double a)
{
 // Purpose: One-step SCAD penalty
  // ----------------------------------------------------------------------
  // Returns lambda * I(x <= l) + I(x > l) * (a*l-x)_+ / (a-1)
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  if(fabs(x) <= l)
    return(l);
  else {
    double x2 = a * l - fabs(x);
  if(x2 < 0)
    return(0);
  return(l * x2 / ((a - 1) * l) );
  }
}

static inline double soft_thr(const double x, const double y) 
{
  // Purpose: Soft thresholding
  // ----------------------------------------------------------------------
  // Returns sign(x)(|x|-y)_+
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen
  
  if(fabs(x) >y)
    {
      if(x>0.0)
	return x - y;
      else
	return x + y;
    }
  return 0;
}

void scale (double *X, const double *weights, const int *n, const int *p,double *mn, double *sd)
{

  // Purpose: Weighted scaling of covariates
  // ----------------------------------------------------------------------
  // Arguments:
  //   X      : data (design matrix flattened by rows)
  //   weights: observation weights - must have sum(weights)=1
  //   n      : number of obs.
  //   p      : number of vars.
  //   mn     : placeholder for mean vector
  //   sd     : placeholder for 'standard deviation'
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  double *msq,*invsd;
  double tmp;
  int npar = *p; 
  int nobs = *n;

  msq =  Calloc(npar, double);
  invsd = Calloc(npar, double);

  for(int i = 0; i < nobs; i++)
    {
      for(int k = 0; k < npar; k++){
	tmp = weights[i] * X[i * npar + k];
	mn[k] += tmp;
	msq[k] += weights[i] * X[i * npar + k] * X[i * npar + k];
      }
    }
  
  for(int k = 0;k < npar; k++){
    tmp = mn[k] * mn[k];
    sd[k] = sqrt(msq[k] - mn[k] * mn[k]);
    if(nobs>1)
      invsd[k]=1/sd[k];
    else
      invsd[k]=1;
  }
  for(int i = 0; i < nobs; i++)
      for(int k = 0; k < npar; k++)
      X[i * npar + k]=(X[i * npar + k] - mn[k]) * invsd[k];

  Free(msq); Free(invsd);
}


void ah (const double *X, const double *inout, const double *tdiff,
	 const double *iatrisk, const int *deathyn, const int *n,const int *p, 
	 double *D, double *d,double *B, const int *univar, const int *usethis,
	 const int *nuse)
{

  // Purpose: Additive hazards regression
  // ----------------------------------------------------------------------
  // Arguments:
  //   X      : data (design matrix flattened by rows)
  //   inout  : indicator, entry (+obs. weight)/exit (-obs. weight)
  //   tdiff  : difference between event times
  //   iatrisk: 1/atrisk
  //   deathyn: is event time a death time?
  //   n      : number of event times
  //   p      : number of vars.
  //   D      : placeholder for D
  //   d      : --"-- d
  //   B      : --"-- B
  //   univar : do univariate analyses?
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  double *H,*Zbar,*covsum,*covsq;
  
  unsigned int npar = *p; 
  
  H = Calloc(npar,double);
  Zbar =  Calloc(npar,double);
  covsum =  Calloc(npar,double);
 
  if(*univar){
    covsq =  Calloc(npar, double);
  } else {
    covsq = Calloc(npar * (*nuse),double);
  }
 /// Loop over episodes
    for(int i = 0; i < *n; i++)
      {
	/// Cumulative covariate processes 
	int ii = 0;
	int kk = 0;

	for(int k = 0; k < npar; k++)
	  {
	    /// If univariate analyses, get only 'variances'
	    int lower = (*univar) * k - 1;

	    H[k] = X[i * npar + k];
	    covsum[k] += H[k] * inout[i];
	    

	    Zbar[k] = covsum[k] * iatrisk[i];
	 
	    // Calculate 'covariances' for this feature?
	    if(usethis[k] > 0)
	      {
		/// Death?
		if( deathyn[i] )
		  {
		    d[k] += inout[i] * (H[k] - Zbar[k]);
		    for(int jj = k; jj > lower; jj--)
		      {
			B[kk] += inout[i] * (H[k] - Zbar[k]) * (H[jj] - Zbar[jj]);
			kk++;
		      }
		  } 
		double tmp=H[k]*inout[i];
		for(int jj = k;jj > lower; jj--)
		  {
		    covsq[ii] += tmp * H[jj];
		    D[ii] += (covsq[ii] - covsum[k] * Zbar[jj]) * tdiff[i];
		    ii++;
		  }
	      }
	  }
      }
    Free(H);
    Free(Zbar);
    Free(covsum);
    Free(covsq);   
}

void ahd (const double *X, const double *inout, const double *tdiff,
	  const double *iatrisk, const int *deathyn, const int *n,const int *p, 
	  double *D, double *d)
{

  // Purpose: Get d and D (~ slimmer version of 'ah')
  // ----------------------------------------------------------------------
  // Arguments:
  //   X      : data (design matrix flattened by rows)
  //   inout  : indicator, entry (+obs. weight)/exit (-obs. weight)
  //   tdiff  : difference between event times
  //   iatrisk: 1/atrisk
  //   deathyn: is event time a death time?
  //   n      : number of event times
  //   p      : number of vars.
  //   D      : placeholder for D
  //   d      : --"-- d
  // ----------------------------------------------------------------------
  // Author: Anders Gorst-Rasmussen

  double H,Zbar;
  double *covsum,*covsq;
  unsigned int npar = *p; 

  covsum = Calloc(npar,double);
  covsq = Calloc(npar,double);

  // Loop over episodes
    for(int i = 0; i < *n; i++)
      {
	for(int k = 0; k < npar; k++)
	  {
	    H = X[i * npar + k];
	    covsum[k] += H * inout[i];

	    Zbar = covsum[k] * iatrisk[i];
	    // Death?
	    if( deathyn[i] )
	      d[k] += inout[i] * (H - Zbar);
	    covsq[k] += H * H * inout[i];
	    D[k] += (covsq[k] - covsum[k] * Zbar) * tdiff[i]; 
	  }
      } 
    Free(covsum);  Free(covsq);
}


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
	  zbar[i * (*p) + k] = covsum[k] * iatrisk[i];
	  bresl[i] -= zbar[i * (*p)+k] * tdiff[i] * beta[k];
	}      
    }  
  Free(covsum);
}


void ahresid (double *start, double *end, double *status, double *X, 
	      double *Zbar, double *times, double *tdiff, double *breslow, 
	      double *beta, int *ntimes,int *p, int *nobs, double *resid)
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
		    resid[i*(*p)+k] += (X[i*(*p)+k] - Zbar[t*(*p)+k]);
		}
	      tmp = 0;
	      for(int k = 0; k < *p; k++)
		tmp += X[i * (*p) + k] * beta[k];
	      for(int k = 0; k < *p; k++)
		resid[i * (*p)+k] += (Zbar[t * (*p)  +k] - X[i* (*p) + k]) * (breslow[t] + tmp * tdiff[t]);
	    }
	}
    }
}

void ahpen (const double * X, const double * inout, const double * wgt, const double * tdiff,
	    const double * times,const double * iatrisk, int * deathyn, const int  *n,const int *p, 
	    double *lambda, const int *nlambda, const double *thresh, const double *penalty,const int *maxit, double *estims, 
	    const int *dfmax,const int *pmax, const double * d, const double * diagD,
	    int *lambdaflag,const double *alpha, int *iter, double *initsol, double *a,
	    const int *nsteps, const double *prefactor)
{

  // Purpose: Penalized estimation for ahaz via CCD
  // -----------------------------------------------
  // Author: Anders Gorst-Rasmussen 




  int npar = *p;   // Number of obs.
  int nobs = *n;   // Number of pars

 
  
  int *check = Calloc(npar, int);             // Check: has variable been considered for inclusion? 
  int *activeidx =  Calloc(npar, int);        // New indices of active variables
  int *activeset =  Calloc(npar,int);         // 0/1-vector: is parameter in current active set? 
  int *oldactiveset =  Calloc(npar, int);     // 0/1-vector: was parameter in last active set? 
  double * betaactive =  Calloc(npar, double);
  double * beta =  Calloc(npar, double);
  double ** D = Calloc(npar,double*);   // Container for




  // Use all variables for active set in first loop
  for(int iii = 0; iii < npar; iii++)
    activeset[iii] = 1;

  int nactive = 0;    // Number of variables which have been considered for inclusion
  double maxdiff,maxtries;
  int test, cnt;

  // Loop over lambda
  for(int l = 0; l < *nlambda; l++)
    {
  
      double lamalph1 = lambda[l] * (*alpha);
      double lamalph2 = (lambda[l] * (1 - *alpha));

      // For multistep penalties
      for(int ll = 0; ll < *nsteps;ll++)
	{
	  maxtries = 50; // Don't restart active sets too many times. Only an issue with exact collinearity.
	  // If activeset != oldactiveset, restart
	  do { 
	    test = 0; // 0/1 - has the active set changed?
	    cnt = 0;   // count - number of convergent passes with current lambda 

	    // One extra loop on all features to obtain "new" active set
	    do { 
	      int i = 0;
	      // Loop over active set until convergence or i>maxit 
	      do {
		maxdiff = *thresh;
		for(int j = 0; j < npar; j++)
		  {
		    if( activeset[j] && penalty[j] > -1 ) 
		      {
			double sum = 0; 
			for(int k = 0; k < nactive; k++)
			  sum += betaactive[k] * D[k][j];
			
			// ------- CCD --------------------
			// Implement penalty of choice here
			double pn = penalty[j] * scad(*prefactor * initsol[j], lamalph1, *a);
			double betanew = soft_thr(d[j] + beta[j] * diagD[j] - sum, pn) /  (diagD[j] + lamalph2);
			// --------------------------------

			if(betanew==0)
			  {
			    // If zero, assign this value + and delete from active set, if ever included
			    beta[j] = 0;
			    activeset[j] = 0;
			    if (check[j]) 
			      betaactive[activeidx[j]] = 0;
			  } else // Include jth feature? 
			    {
			      if(check[j] != 0) // First time betanew is included? If yes, get row in D
				betaactive[activeidx[j]] = betanew; 
			      else {
				  betaactive[nactive] = betanew;
				  activeidx[j] = nactive;   // New index of jth feature
				  check[j] = 1;
				  
				  // Calculate row in D
				  // Messy - but, we want few calculations in innermost loop
				  D[nactive] =  Calloc(npar, double); 
				  double xx[nobs][2];
				  double covsum = 0;
				  for(int ii = 0; ii < nobs; ii++)
				    {
				      covsum += X[j * nobs + ii] * inout[ii];
				      xx[ii][1] = covsum * iatrisk[ii] * tdiff[ii];
				      xx[ii][2] = X[j * nobs + ii]* times[ii] * wgt[ii];
				    }	     
				  for(int jj = 0; jj < npar; jj++)
				    {
				      double sum1 = 0; double sum2 = 0;
				      for(int ii = 0; ii < nobs; ii++)
					{
					  sum1 += X[jj * nobs + ii] * inout[ii];
					  sum2 += X[jj * nobs + ii] * xx[ii][2] - xx[ii][1] * sum1;
					}
				      D[nactive][jj] = sum2;
				    }
				  // If too many active features, break out
				  if( nactive >= *pmax ){
				    *lambdaflag = l;
				    goto out;  
				  }
				  nactive++;
			      }
		     
			      // Change from last beta - convergence?
			      double bdiff = fabs((beta[j] - betanew) /  beta[j]);
			      if(bdiff > maxdiff)
				maxdiff = bdiff;	
			      			      
			      // Assign updated parameter
			      beta[j] = betanew;
			    }  
	
		      }
		  }
		i++;
	      }  while(maxdiff> *thresh && i < *maxit);
	      
	      // Keep track of max number/iterations
	      if (i >= *iter)
		*iter = i;
	      
	      // If we have done one iteration, we need one more complete iteration
	      if (cnt == 0)
		{
		  for(int jj = 0; jj < npar; jj++){
		    oldactiveset[jj] = activeset[jj];
		    activeset[jj] = 1;
		  }
		}  
	      cnt++;} while(cnt<2);
	    
	    // Has active set changed - if not, start over.
	    for(int jj = 0; jj < npar; jj++)
	      if(oldactiveset[jj] != activeset[jj])
		{
		  test++;
		  maxtries++;
		  break;
		}
	  } while(test!=0 && maxtries);
	  
	  // Update 'initial solution' if multistep algorithm
	  if(*nsteps > 1)
	    for(int jj = 0; jj < npar; jj++)
	      initsol[jj] = beta[jj];
	}
      // Save relevant data to placeholder
      int params = 0;
      for(int jj = 0; jj < npar; jj++){
	params += (beta[jj] != 0);
	estims[l * npar + jj] = beta[jj];
      }
      // Tell R at which lambda to stop
      if(params > *dfmax)
	{
	  *lambdaflag = l+1;
	  break;
	}
    }
  
 out:;
  
  // Cleanup
  for(int j = 0; j < nactive; j++)
    Free(D[j]);
  Free(D);
  Free(check);Free(activeidx);
  Free(activeset); Free(oldactiveset);
  Free(beta);Free(betaactive);
}
