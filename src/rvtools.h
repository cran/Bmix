/*  A set of basic random variable tools that work with 
 *  the randomkit.h Mersenne Twister pseudo-RNG.
 *  
 */

#ifndef __RVTOOLS_H__
#define __RVTOOLS_H__

#include <stdio.h>
 
void* newRNGstate(unsigned long s);
void* newRNGstate_rand(void *s);
void deleteRNGstate(void *seed);
unsigned long three2lstate(int *state);

double runi(void *state); 
double rnor(void *state);
double rst(double nu, double m, double s2, void *state);
double rexpo(double lambda, void *state);
double rgamma1(double aa, void *state);
double rgamma2(double alpha, void *state);
double rgam(double alpha, double beta, void *state);
void rdir(double *theta, double *eta, int d, void *state);
double rbet(double aa, double bb, void *state);

int indexsample(int *ind, int n, int num_probs, double *probs, void *state);
int indexdraw(int num_probs, double *probs, void *state);
void worsample(int *ind, int n, int num_probs, double *probs, void *state);

double ldNbb(double y, double a, double b, double r);
double ldPg(double y, double a, double b, double r);
double ldbinbet(double y, double a, double b, double n);
double ldst(double z, double nu, double m, double s2);

double pbar(double v, double a, double b, double p, void*state);
double stickbreak(double *w, int L, double sticka, void *state);

double logit(double x);
double expit(double z);

#endif
