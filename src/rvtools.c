#include <math.h>
#include <stdlib.h>
#include <assert.h> 
#include <Rmath.h>

#include "rvtools.h"
#include "latools.h"
#include "randomkit.h"

/* 
 * newRNGstate:
 * 
 * seeding the random number generator,
 */

void* newRNGstate(unsigned long s)
{
   rk_state* state = (rk_state*) malloc(sizeof(rk_state));
   rk_seed(s, state);
   return (void*) state;
}


/*
 * newRNGstate_rand:
 *
 * randomly generate a new RNG state based on a random draw from the
 * current state
 */

void* newRNGstate_rand(void *s)
{
  unsigned long lstate;
  int state[3];
  state[0] = 100*runi(s);
  state[1] = 100*runi(s);
  state[2] = 100*runi(s);
  lstate = three2lstate(state);
  return(newRNGstate(lstate));
}


/*
 * three2lstate:
 *
 * given three integers (positive) , turning it into
 * a long-state for the RNG seed
 */

unsigned long three2lstate(int *state)
{
  unsigned long lstate;
  assert(state[0] >= 0);
  assert(state[1] >= 0);
  assert(state[2] >= 0);
  lstate = state[0] * 1000000 + state[1] * 1000 + state[2];
  return(lstate);
}

      
/*
 * deleteRNGstate:
 *
 * free memory for RNG seed
 */

void deleteRNGstate(void *state)
{
   free((rk_state*) state);
}




/* 
 * runi:
 * 
 * one from a uniform(0,1)
 */

double runi(void *state)
{
    unsigned long rv;
    assert(state);
    rv = rk_random((rk_state*) state);
    return ((double) rv) / RK_MAX;
}


/*
 * rnor:
 * 
 * one draw from a standard normal 
 */

double rnor(void *state)
{
  double z;
  z = rk_gauss((rk_state*) state);
  return z;
}

/*
 * rnor:
 * 
 * one draw from a student's t
 */

double rst(double nu, double m, double s2, void *state)
{

  double z = rnor(state);
  double x = rgam(0.5*nu, 0.5, state);
  return m + z*sqrt(s2*nu/x);

}

double rexpo(double lambda, void *state)
/*
 * Generates from an exponential distribution
 */
{
    double random, uniform;
    uniform = runi(state);
    random = 0.0 - (1/lambda) * log(uniform);
    return random;
}



/*
 * rgamma1:
 * 
 * Generates a draw from a gamma distribution with alpha < 1
 * 
 */

double rgamma1(double alpha, void *state)
{
  double uniform0, uniform1;
  double random, x;
  
  /* sanity check */
  assert(alpha > 0);
  
  /* int done = 0; */
  uniform0 = runi(state);
  uniform1 = runi(state);
  if (uniform0 > M_E/(alpha + M_E))
    {
      random = 0.0 -log((alpha + M_E)*(1-uniform0)/(alpha*M_E));
      if ( uniform1 > pow(random,alpha - 1))
	return -1;
      else 
	return random;
    }
  else
    {
      x = (alpha + M_E) * uniform0 / M_E;
      random = pow(x,1/alpha);
      if ( uniform1 > exp(-random))
	return -1;
      else
	return random;
    } 
}


/*
 * rgamma2:
 * 
 * Generates a draw from a gamma distribution with alpha > 1
 *
 * 
 */

double rgamma2(double alpha, void *state)
{
  double uniform1,uniform2;
  double c1,c2,c3,c4,c5,w;
  double random;
  int done = 1;
  
  /* sanity check */
  assert(alpha > 0);
  
  c1 = alpha - 1;
  c2 = (alpha - 1/(6 * alpha))/c1;
  c3 = 2 / c1;
  c4 = c3 + 2;
  c5 = 1 / sqrt(alpha);
  do
    {
      uniform1 = runi(state);
      uniform2 = runi(state);
      if (alpha > 2.5)
        {
	  uniform1 = uniform2 + c5 * (1 - 1.86 * uniform1);
        }
    }
  while ((uniform1 >= 1) || (uniform1 <= 0));
  
  w = c2 * uniform2 / uniform1;
  if ((c3 * uniform1 + w + 1/w) > c4)
    {
      if ((c3 * log(uniform1) - log(w) + w) >= 1)
        {
	  done = 0;
        }
    }
  if (done == 0)
    return -1;
  random = c1 * w; 
  return random;
}


/*
 * rgam:
 * 
 * Generates from a general gamma(alpha,beta) distribution
 * Parametrization as in the Gelman's book ( E(x) = alpha/beta )
 */

double rgam(double alpha, double beta, void *state)
{
  double random = 0;
  
  /* sanity checks */

  assert(alpha>0 && beta>0);
  
  if (alpha < 1)
    do {
      random = rgamma1(alpha, state)/beta; 
    } while (random < 0 );
  if (alpha == 1)
    random = rexpo(1.0, state)/beta; 
  if (alpha > 1)
    do {
      random = rgamma2(alpha, state)/beta; 
    } while (random < 0);
  return random;
}


/*
 * rdirichlet:
 *
 * Generates from a dirichlet distribution.
 * 
 */


void rdir(double *theta, double *eta, int d, void *state)
{
  double *x = new_dvec(d);
  unsigned int i;
  double xsum = 0.0;

  for(i=0; i<d; i++){ 
    assert(eta[i] > 0);
    x[i] = rgam(eta[i], 1.0, state);
    xsum += x[i];
  }
  for(i=0; i<d; i++) theta[i] = x[i]/xsum;

  free(x);
}


/*
 * rbet:
 * 
 * one random draw from the beta distribution
 * with parameters alpha and beta.
 */

double rbet(alpha, beta, state)
double alpha, beta;
void *state;
{
   double g1,g2;
   assert(alpha>0);
   assert(beta>0);
   g1 = rgam(alpha, 1.0, state);
   g2 = rgam(beta, 1.0, state);

   return g1/(g1+g2);
}

/*
 * Without-replacement sampling
 *
 */
 
void worsample(int *ind, int n, int num_probs, double *probs, void *state){
  int i, j;
  double *dprobs = new_dup_dvec(probs, num_probs);
  int *index = new_iseq(0, num_probs-1);
  int dnp = num_probs;
  double pi;
  int tempi;
  for(i=0; i<n; i++){
    indexsample(&tempi, 1, dnp, dprobs, state);
    ind[i] = index[tempi];
    pi = dprobs[tempi];
    dnp -=1;
    for(j=0; j<tempi; j++){ dprobs[j] *= 1.0/(1.0-pi); }
    for(j=tempi; j<dnp; j++){ 
      dprobs[j] = dprobs[j+1]/(1.0-pi); 
      index[j] = index[j+1];
    }
  }

  free(index);
  free(dprobs);
}


/* 
 * indexsample:
 *
 * Just simple "with replacement" index sampling
 * Returns the number of unique components.
 * 
 */

int indexsample(int *ind, int n, int num_probs, double *probs, void *state)
{
  double pick;
  int i, counter;
  double *cumprob = new_dvec(num_probs);
  double *selected = new_dzero(num_probs);
  
  assert(num_probs > 0);
  assert(n > 0);
  
  assert(probs[0] >= 0);
  cumprob[0] = probs[0];
  for(i=1; i<num_probs; i++) {
    assert(probs[i] >= 0);
    cumprob[i] = cumprob[i-1] + probs[i];
  }
  if(cumprob[num_probs-1] < 1.0) cumprob[num_probs-1] = 1.0;
  
  for(i=0; i<n; i++) {
    counter = 0;
    pick=runi(state);
    while(cumprob[counter] < pick) counter++;
    ind[i] = counter;
    selected[counter]++;
  }

  counter = 0;
  for(i=0; i<num_probs; i++) if(selected[i] > 0) counter++; 
  free(cumprob);
  free(selected);
  return counter;
}

int indexdraw(int num_probs, double *probs, void *state){
  int ind;
  indexsample(&ind, 1, num_probs, probs, state);
  return ind;
}


double ldNbb(double y, double a, double b, double r)
{
  double numerator = lgammafn(a+b) + lgammafn(a+r) + lgammafn(b+y);
  double denominator = lgammafn(a) + lgammafn(b) + lgammafn(a+b+r+y);
  return numerator - denominator + lchoose(r+y-1, r-1);
}


double ldPg(double y, double a, double b, double r)
{
  double ld = lgammafn(a+y) - lgammafn(a) - lgammafn(y+1.0);
  ld += y*log(r) + a*log(b) - (a+y)*log(b+r);
  return ld;
}


double ldbinbet(double y, double a, double b, double n)
{
  double numerator = lgammafn(a+b) + lgammafn(a+y) + lgammafn(b+n-y);
  double denominator = lgammafn(a) + lgammafn(b) + lgammafn(a+b+n);
  return numerator - denominator + lchoose(n, y);
}


double ldst(double z, double nu, double m, double s2)
{
  double ld = lgammafn((nu+1.0)*0.5) - lgammafn(nu*0.5);
  ld +=  -log(sqrt(nu*s2)) - M_LN_SQRT_PI;
  ld += -0.5*(nu+1.0)*log( 1.0 + (z-m)*(z-m)/(nu*s2) );
  return ld;
}


double pbar(double v, double a, double b, double p, void*state)
{
  assert(a>p);
  double u = rbet(b,a-p, state);
  double w = rbet(p,a-p, state);
  return 1.0 - u*(1-w*v);
}

double stickbreak(double *w, int L, double sticka, void *state)
{  
  double *zeta = new_dvec(L);
  int l;

  for(l=0; l<L; l++) zeta[l] = rbet(1.0, sticka, state);

  w[0] = zeta[0];
  double wsum = w[0];

  for(l=0; l<(L-1); l++){
    w[l+1] = ((1-zeta[l])/zeta[l])*w[l]*zeta[l+1];
    wsum += w[l+1];
  }

  for(l=0; l<L; l++) w[l] = w[l]/wsum; 

  free(zeta);
  return 1.0-wsum;
}


double logit(double x)
{
  if(x <= 0.0 || x >= 1.0){ printf("bad x in logit\n"); return 0.0; }
  else return log(x) - log(1-x);
}

double expit(double z)
{
  return exp(-log(1+exp(-z)));
} 
