/*  Written by Matt Taddy, Chicago Booth  */

#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <vector>
#include "matrix.h"
using namespace std;


class Particle
{
 private:
  /* dims and nums */
  int time;
  int obs;
  int dim;
  double d;
  int m;

  vector<double> n;
  vector<int> k;
  vector<double> v;
  Matrix p;

  /* params  */
  double alpha;
  double rho;
  Matrix lambda;
  double kappa;
  double nu;
  Matrix Omega;
    
  /* categorical items */
  int cat;
  vector<int> levels;
  vector<Matrix> aQ; 
  vector< vector<Matrix> > R;

  /* hyperparams */
  double gamO;
  Matrix psiO;

  /* data */
  vector<Matrix> ZM;

  /* suff-stats */
  vector<Matrix> zbar;
  vector<Matrix> S;

  /* post pred params */
  // for each component: 
  vector<Matrix> a;
  vector<Matrix> B;
  vector<double> c;
  vector<Matrix> D;
  // and for the base 
  Matrix a0;
  Matrix B0;
  double c0;

  /* stats for augment step */
  vector<double> ct;
  int mbk;
  vector<double> vbk;
  double pct;
  

 public:

  Particle(int *dims, double *params);
  /* NO POINTERS <-> use synthesized copy control */
  // Particle(const Particle &p);
  // ~Particle(void);
  // Particle& operator=(const Particle &p);

  void Params(double *params);
  void ABCD(int j);

  void Push(double *zm);
  void Erase(int j);
  void Add(int i);
  void Remove(int i);
  void CalcP();

  Matrix Probs(double *zm);
  double PostPred(double *zm);
  double Augment(int newtime, void *state);
  void Propagate(double *zm,  void *state);
  void DrawG0(void *state);
  void DrawK(int i, void *state);
  void DrawFull(void *state); // includes the above, but not V

  int sumN();
  int sumCnt();
  int getM();
  void writeK(int *out);

  void Print(int ind);
  void Print(int ind, int tt);
  void Read(int ind, int num);

};

#endif
  

 
