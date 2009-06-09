/*  Written by Matt Taddy, Chicago Booth  */

extern "C" {
#include "rvtools.h"
#include "latools.h"
#include <stdio.h>
}
#include "particle.h"
#include <assert.h>
#include <Rmath.h>

using namespace std;

Particle::Particle(int *dims, double *params)
{
  time = 0;
  obs = 0;
  pct = 1.0;
  dim = dims[0];
  d = (double) dim;
  cat = dims[1];
  if(cat > 0) for(int r=0; r<cat; r++) levels.push_back(dims[2+r]); 
    
  Params(params);
  m = 0;
  ABCD(m); 
  CalcP();
}

void Particle::Params(double *params)
{
  alpha = params[0]; assert(alpha > 0.0);
  rho = params[1]; assert(rho >= 0 && rho <= 1);
  lambda = Matrix(dim, 1, &params[2]); 
  kappa = params[dim+2]; assert(kappa > 0.0);
  nu = params[dim+3]; assert(nu > 0.0);
  gamO = params[dim+4];
  psiO = Matrix(dim, dim, &params[dim+5], true);
  Omega = psiO*gamO;
  for(int i=0; i<cat; i++)
    aQ.push_back(Matrix(levels[i], 1, 1.0/((double) levels[i])));
}

void Particle::ABCD(int j)
{ 
  if(j==m){ a0 = lambda;
    c0 =  2.0*nu - d + 1.0;
    B0 = Omega*( 2.0*(kappa+1.0)/(kappa*c0) ); 
    return; }
  a[j] = ( kappa*lambda + n[j]*zbar[j] )/( kappa + n[j] );
  c[j] = 2.0*nu + n[j] - d + 1.0;
  D[j] =  S[j] + (lambda-zbar[j])*Transpose(lambda-zbar[j])*kappa*n[j]/( kappa + n[j] ); 
  B[j] = 2.0*( kappa + n[j] + 1.0 )*(Omega + 0.5*D[j])/( ( kappa + n[j] )*c[j] );
  B[j].SetSym(true); 
  D[j].SetSym(true);
}

void Particle::Push(double *zm)
{
  n.push_back(1);
  zbar.push_back(Matrix(dim, 1, zm));
  if(v.size() > 0)
    { v.push_back( rbeta( 2.0, alpha ) ); // uses R's rng
      ct.push_back( 1.0 ); }
  if(cat>0) 
    { vector<Matrix> catmat;
      for(int r=0; r<cat; r++)
	{ catmat.push_back(Matrix(levels[r],1));
	   catmat[r][0][(int) zm[dim+r] ] = 1.0;  } 
      R.push_back(catmat); }
  S.push_back(Matrix(dim,dim));
  a.push_back(Matrix(dim, 1));
  B.push_back(Matrix(dim,dim));
  c.push_back(0.0);
  D.push_back(Matrix(dim,dim));
  m++;
  CalcP();
  ABCD(m-1); 
}

void Particle::Erase(int j){
  if(v.size()>0)
    { cout << "Shouldn't be running MCMC with ddp weights; reset." << endl;
      v.clear(); ct.clear(); }
  m+= -1;
  assert(n.size() == m+1 && n[j]==0);
  n.erase(n.begin()+j);  
  CalcP();
  zbar.erase(zbar.begin()+j);
  S.erase(S.begin()+j);
  if(cat>0) R.erase(R.begin()+j);
  a.erase(a.begin()+j);
  B.erase(B.begin()+j);
  c.erase(c.begin()+j);
  D.erase(D.begin()+j);
  for(int i=0; i<obs; i++){ assert(k[i] != j); if(k[i]>j) k[i]--; }
}

void Particle::Add(int i){   
  if(k[i] == m) Push(ZM[i][0]);
  else{
    Matrix zi = Matrix(dim, 1, ZM[i][0]);
    int j = k[i];
    n[j]++;

    if(v.size()==0) CalcP(); 
    else ct[j]++;

    S[j] += zi*Transpose(zi) + (n[j]-1.0)*zbar[j]*Transpose(zbar[j]);
    zbar[j] = ( (n[j]-1.0)*zbar[j] + zi )/n[j];
    S[j] += -n[j]*zbar[j]*Transpose(zbar[j]); 
    ABCD(j);
    for(int r=0; r<cat; r++) R[j][r][0][(int) ZM[i][0][dim+r] ] += 1.0;  
  }
}

void Particle::Remove(int i){
  if(v.size()>0)
    { cout << "Shouldn't be running MCMC with ddp weights; reset." << endl;
      v.clear(); }
  int j = k[i];
  k[i] = -1;
  n[j]+= -1.0;
  if(n[j] == 0.0) Erase(j);
  else{
    Matrix zi = Matrix(dim, 1, ZM[i][0]);
    CalcP();
    S[j] +=  (n[j]+1.0)*zbar[j]*Transpose(zbar[j])-zi*Transpose(zi);
    zbar[j] = ( (n[j]+1.0)*zbar[j] - zi )/n[j];
    S[j] += -n[j]*zbar[j]*Transpose(zbar[j]); 
    ABCD(j);
    for(int r=0; r<cat; r++) R[j][r][0][(int) ZM[i][0][dim+r] ] -= 1.0;
  }
}

void Particle::CalcP(void){ 
  double sumv;
  p = Matrix(m+1, 1);
  if(v.size()==0){
    sumv = alpha;
    for(int j=0; j<m; j++) sumv += n[j];
    for(int j=0; j<m; j++) p[0][j] = n[j]/sumv;
    p[0][m] = alpha/sumv; 
  }
  else{
    sumv = 1.0;
    for(int j=0; j<m; j++)
      { p[0][j] = sumv*v[j]; sumv -= p[0][j]; }
    p[0][m] = sumv;
  }
}

Matrix Particle::Probs(double *zm)
{
  const Matrix ZZ = Matrix(dim, 1, zm);
  Matrix pzm = Matrix(m+1, 1);
  pzm[0][m] = exp( log(p[0][m]) + ZZ.ldST(c0, a0, B0) );
  for(int j=0; j<m; j++) 
    pzm[0][j] = exp( log(p[0][j]) + ZZ.ldST(c[j], a[j], B[j]) );
  if(cat > 0) 
    for(int r=0; r<cat; r++){
      pzm[0][m] = exp(log(pzm[0][m]) + log(aQ[r][0][(int) zm[dim+r]]) - log(aQ[r].Sum()) );
      for(int j=0; j<m; j++) 
	pzm[0][j] = exp(log(pzm[0][j]) 
			+ log( aQ[r][0][(int) zm[dim+r]] + R[j][r][0][(int) zm[dim+r]] )
			- log( aQ[r].Sum() + R[j][r].Sum() ) );
        }

  return pzm;
}

double Particle::PostPred(double *zm){ Matrix pzm = Probs(zm); return pzm.Sum(); }


void Particle::Propagate(double *zm, void *state)
{
  ZM.push_back(Matrix(dim+cat,1,zm));
  Matrix pk = Probs(zm);
  pk.Normalize();
  // sample k and propagate
  k.push_back(indexdraw(m+1, pk[0], state));
  Add(obs);
  obs++;   
}

void Particle::DrawG0(void *state)
{
  /* draw latent kenel params */
  vector<Matrix> Sigmai;
  for(int j=0; j<m; j++){
    Sigmai.push_back(Matrix(dim, dim));
    Matrix var = Omega + D[j];
    Sigmai[j].rWSH( (int) (nu + n[j]), var, state);
  }

  /* omega */
  Matrix SS = Inverse(psiO);
  for(int j=0; j<m; j++) SS += Sigmai[j];
  double df = nu*((double) m) + gamO;
  Omega.rWSH( ((int) df), SS, state );
  /* re-calculate postpred params */
  for(int j=m; j>=0; j--) ABCD(j);
}


void Particle::DrawK(int i, void *state)
{   
  Remove(i);  // take Zi out

  // build the probabilities
  Matrix pk = Probs(ZM[i][0]);
  pk.Normalize();

  // propose k
  indexsample(&k[i], 1, m+1, pk[0], state);

  Add(i);  // put Zi back
}

void Particle::DrawFull(void *state)
{
  for(int i=0; i<obs; i++) DrawK(i, state);
  DrawG0(state);
}

double Particle::Augment(int newtime, void *state){
  /* draw initial probs */ 
  if(v.size()==0)
      { mbk = m;
	double nsum = 0.0; 
	for(int j=0; j<m; j++) nsum += n[j];
	for(int j=0; j<m; j++)
	  { nsum -= n[j];
	    vbk.push_back(rbet(1.0+n[j], alpha+nsum, state));
	    ct.push_back(n[j]); } 
	v = vbk; }

  /* newtime! */
  
  if(time != newtime){ 
    assert(time+1 == newtime);
    time = newtime; 
    mbk = m;
    vbk = v; 
    for(int j=0; j<m; j++) ct[j] = 0.0;  }
  
  /* propate from the pbar prior for j<mbk*/
  pct = 1.0;
  double *sct = new_dvec(m);
  sct[m-1] = ct[m-1]; 
  for(int j=m-2; j>=0; j--) sct[j] = sct[j+1]+ct[j];

  for(int j=0; j<mbk; j++)
    { v[j] = pbar(vbk[j], 1.0, alpha, rho, state); 
      pct = exp(log(pct) + dbinom(ct[j], sct[j], v[j], 1)); }
  for(int j=mbk; j<m; j++)
    { v[j] = rbet(1.0+ct[j], alpha+sct[j]-ct[j], state);  }

  if(sct[0] == 0.0) assert(pct = 1.0);
  free(sct);
  CalcP(); 
  return pct;
}

int Particle::sumN(void)
{ double sn =0.0;  
  for(int j=0; j<m; j++) sn += n[j]; 
  return (int) sn; }

int Particle::getM(void){ return m; }

void Particle::writeK(int *out){ for(int i=0; i<k.size(); i++) out[i] = k[i];  }

void Particle::Print(int ind, int tt){
  int temp = time;
  time = tt;
  Print(ind);
  time = temp;
}

void Particle::Print(int ind)
{
 /* print the results to file*/
  FILE *file;
  char format[] = ".particle%d.%d.%g.txt";
  char filename[sizeof format+100];
  sprintf(filename,format,ind,time+1,rho);
  file = fopen(filename,"w"); 
  assert(file);
      
  /* print the base measure and prior pred */
  fprintf(file, "%g ", alpha);
    for(int i=0; i<dim; i++) fprintf(file, "%g ", lambda[0][i]);
    for(int j=0; j<dim; j++) 
      for(int i=0; i<dim; i++) 
	fprintf(file, "%g ", Omega[j][i]);
    for(int r=0; r<cat; r++)
      for(int s=0; s<levels[r]; s++)
	fprintf(file, "%g ", aQ[r][0][s]);
    fprintf(file, "%g ", p[0][m] );
    for(int i=0; i<dim; i++) fprintf(file, "%g ", a0[0][i]);
    for(int j=0; j<dim; j++) 
      for(int i=0; i<dim; i++) 
	fprintf(file, "%g ", B0[j][i]);
    fprintf(file, "%g ", c0);
    fprintf(file, "\n");

  /* print suff stat and momments for each component */
  for(int l=0; l<m; l++){
    fprintf(file, "%g ", n[l]);
    for(int i=0; i<dim; i++) fprintf(file, "%g ", zbar[l][0][i]);
    for(int j=0; j<dim; j++) 
      for(int i=0; i<dim; i++) 
	fprintf(file, "%g ", S[l][j][i]);
    for(int r=0; r<cat; r++)
      for(int s=0; s<levels[r]; s++)
	fprintf(file, "%g ", R[l][r][0][s]);
    fprintf(file, "%g ", p[0][l] );
    for(int i=0; i<dim; i++) fprintf(file, "%g ", a[l][0][i]);
    for(int j=0; j<dim; j++) 
      for(int i=0; i<dim; i++) 
	fprintf(file, "%g ", B[l][j][i]);
    fprintf(file, "%g ", c[l]);
    fprintf(file, "\n");
  }
  fclose(file);
}

void Particle::Read(int ind, int num)
{
/* print the results to file*/
  char format[] = ".particle%d.%d.txt";
  char filename[sizeof format+100];
  sprintf(filename,format,ind,num);
  FILE *file = fopen(filename,"r");
  if(!file)
    { cout << "Missing file '.particle" << ind << "." << num << ".txt' for input." << endl;
      return; }

  while(fgetc(file) != '\n'){ /* skip first line of prior params */ }

  /* read in the compnents */
  double nj; 
  int suffdim = dim*(dim+1);
  for(int r=0; r<cat; r++) suffdim += levels[r];
  double *suffstats = new_dvec(suffdim);
  while(fscanf(file, "%lf", &nj)==1){ 
    n.push_back(nj);
    for(int i=0; i<suffdim; i++)
      if(fscanf(file, "%lf", &suffstats[i]) != 1) cout << "read error" << endl;
    /* now, push back the suffstats */
    int pnt0=0;
    zbar.push_back(Matrix(dim, 1, suffstats));
    pnt0 += dim;
    S.push_back(Matrix(dim,dim, &suffstats[pnt0]));
    pnt0 += dim*dim;
    if(cat>0) 
      { vector<Matrix> catmat; 
	for(int i=0; i<cat; i++)
	  { catmat.push_back(Matrix(levels[i],1, &suffstats[pnt0]));
	    pnt0 += levels[i]; } 
	R.push_back(catmat); }
    a.push_back(Matrix(dim, 1));
    B.push_back(Matrix(dim,dim));
    c.push_back(0.0);
    D.push_back(Matrix(dim,dim));
    m++;
    ABCD(m-1);
    /******** end push **********/  
    while( fgetc(file) != '\n'){  /* just read to the end of line */ }
  } 
  free(suffstats);
  
  CalcP(); // v is empty

  fclose(file);
}


