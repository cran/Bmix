/*  Written by Matt Taddy, Chicago Booth  */

#include <Rmath.h>
#include <R.h>
extern "C" {
#include "rvtools.h"
#include "latools.h"
#include "rhelp.h"
}
#include "matrix.h"
#include "particle.h"
#include "assert.h"

#include <stdio.h>
#include <iomanip>
#include <vector>
#include <math.h>

using namespace std;

extern "C" {

void mixsample(int* rstate,
	       int* iopar, // read?, print?
	       int* nums,  // n.particle, n.iterations
	       int* dims,  // dim, cat, levels
	       int* Total,
	       int* time,
	       double* Z,
	       double* params, // interpreted in particle
	       double* margllhdout,
	       double* m_out,
	       int* k_out
	       ) 
{  

  double rho = params[1];

  /* initiate the twister with a state from R */
  void* state = NULL;
  unsigned int lstate = three2lstate(rstate);
  state = newRNGstate(lstate);

  /* What's on the menu? */
  int READ = iopar[0];
  int PRINT = iopar[1];
  int N = nums[0];
  int niter = nums[1];
  /* import the data and dimensions */
  int d = dims[0]+dims[1];
  int total = *Total;

  /* build and initiate the particle sets */
  vector<Particle> pset, resampled;
  for(int i=0; i<N; i++) pset.push_back(Particle(dims, params)); 
  if(READ>0) for(int i=0; i<N; i++) pset[i].Read(i+1, READ); 
  

  /* variables and vectors used during filtering */
  int *index = new_ivec(N);
  double *prob = new_dvec(N);
  double *weight = drep(1.0/((double) N), N);
  double var, ncomp;
  int np = 0;

  double pred;
  double mllhd = 0.0;
  double mo;


  /* Filter the observations */
  for(int r=0; r<total; r++)
    { 
	
      /* Augment */
      if(time[r] > 0)
	if(rho < 1.0 && rho > 0.0){
	  for(int i=0; i<N; i++)
	    { if(time[r] != time[r-1] && PRINT) pset[i].Print(i+1);
	      weight[i] = pset[i].Augment(time[r], state); }
	  normalize(weight, N); }// weights sum to one
	else
	  if(rho == 0.0 && time[r] != time[r-1]){
	    if(PRINT) for(int i=0; i<N; i++) pset[i].Print(i+1, time[r-1]);
	      pset.clear();
	      for(int i=0; i<N; i++) pset.push_back(Particle(dims, params)); }

      /* Resample */
      ncomp = 0.0;
      pred = 0.0;
      for(int i=0; i<N; i++)
	{ prob[i] = weight[i]*pset[i].PostPred(&Z[r*d]);  
	  ncomp += (double) pset[i].getM();
	  pred += prob[i]; }
      ncomp = ncomp/((double) N);
      var = normalize(prob, N);
      np = indexsample(index, N, N, prob, state);

      /* track the marginal likelihood */
      mllhd += log(pred);
      margllhdout[r] = mllhd;
       
      for(int i=0; i<N; i++) resampled.push_back(pset[index[i]]);
      pset.clear(); pset.swap(resampled);
      assert( pset.size()==N && resampled.size()==0 );

      /* Propagate */
      for(int i=0; i<N; i++) pset[i].Propagate(&Z[r*d], state);
      for(int i=0; i<N; i++) pset[i].DrawG0(state);  
      
      /* What just happened? */
      if(r %100 == 0)
		bobbys_printf(bobbys_stdout, "r = %d, t = %d, resampled = %d, avg(clusters) = %g\n",r,time[r],np,ncomp); 

      /* Track the number of components */  
      mo = 0.0;
      for(int i=0; i<N; i++) mo += (double) pset[i].getM(); 
      m_out[r] = mo/((double) N); 

    }

  /* Run MCMC and record results */
  
  if(niter >10){
    if(N > 1)
      bobbys_printf(bobbys_stdout, "More than 10 iterations: I'll run MCMC for 1st particle only.");
    for(int j=0; j<niter; j++)
      { pset[0].DrawFull(state);
	pset[0].writeK(&k_out[j*total]);
	if(PRINT) pset[0].Print(j+1); }
  }
  else{
    for(int j=0; j<niter; j++)
      for(int i=0; i<N; i++) if(j>0) pset[i].DrawFull(state); 
    if(k_out) for(int i=0; i<N; i++) pset[i].writeK(&k_out[i*total]); 
    if(PRINT) for(int i=0; i<N; i++) pset[i].Print(i+1);
  }    

    
  /* clean-up */
  free(prob); free(index); free(weight); 
  deleteRNGstate(state); state=NULL;

}

void rwish(int* rstate, int *dim, int *nu, double *psi){

  /* initiate the twister with a state from R */
  void* state = NULL;
  unsigned int lstate = three2lstate(rstate);
  state = newRNGstate(lstate);
	   
  Matrix Psi = Matrix(*dim, *dim, psi, true);
  Matrix Omega = Matrix(*dim, *dim);

  Omega.rWSH(*nu, &Psi, state);
  for(int i=0; i<(*dim); i++)
    copy_dvec(&psi[i*(*dim)], Omega[i], (*dim));

  deleteRNGstate(state); state=NULL;
}
  


}//extern

