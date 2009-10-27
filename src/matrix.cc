extern "C" {
#include "rvtools.h"
#include "latools.h"
}

#include "matrix.h"
#include <Rmath.h>
#include <R.h>
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

Matrix::Matrix(void)
{
  sym=false;
  nrow=0;
  ncol=0;
  M = NULL;
  //cout << "Null matrix call" << endl;
}

Matrix::Matrix(int nrow_in, int ncol_in) // default to a zero matrix
{
  sym = false;
  nrow = nrow_in;
  ncol = ncol_in;
  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = 0.0;
}


Matrix::Matrix(int nrow_in, int ncol_in, double val) // default to a zero matrix
{
  sym = false;
  nrow = nrow_in;
  ncol = ncol_in;
  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = val;
}

Matrix::Matrix(int nrow_in, int ncol_in, const double *vals)
{
  sym = false;
  nrow = nrow_in;
  ncol = ncol_in;

  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = vals[j*nrow+i];
}

Matrix::Matrix(int nrow_in, int ncol_in, const double *vals, bool sym_in)
{
  sym = true;
  nrow = nrow_in;
  ncol = ncol_in;
  if(nrow!=ncol) error("Trying to declare a symmetric matrix with nrow!=ncol.");

  M = new_mat(nrow, ncol);

  for(int j=0; j<ncol; j++){
    for(int i=j; i<nrow; i++){
      assert(vals[j*nrow+i]==vals[i*nrow+j]);
      M[j][i] = M[i][j] = vals[j*nrow+i]; 
    }
  }
}


Matrix::Matrix(int nrow_in, int ncol_in, const Matrix &mat, int *rows, int *cols)
{
  sym = false;
  nrow = nrow_in;
  ncol = ncol_in;

  M = new_mat(nrow, ncol);
 
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = mat.M[cols[j]][rows[i]];
}


Matrix::Matrix(int nrow_in, int ncol_in, const Matrix *mat, int *rows, int *cols)
{
  sym = false;
  nrow = nrow_in;
  ncol = ncol_in;

  M = new_mat(nrow, ncol);
 
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = mat->M[cols[j]][rows[i]];

}


Matrix::Matrix(const Matrix &mat)
{
  sym = mat.sym;
  nrow = mat.nrow;
  ncol = mat.ncol;
  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = mat.M[j][i];
}
  


Matrix::Matrix(const Matrix *mat)
{
  sym = mat->sym;
  nrow = mat->nrow;
  ncol = mat->ncol;
  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = mat->M[j][i];
}

Matrix::~Matrix(void)
{
  delete_mat(M);
}

Matrix& Matrix::operator=(const Matrix &rhs)
{  
  sym = rhs.sym;
  nrow = rhs.nrow;
  ncol = rhs.ncol;
  double **MC = new_dup_mat(nrow, ncol, rhs.M); // to be safe for self assignment
  if(M) delete_mat(M); 
  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = MC[j][i];
  delete_mat(MC);
  return *this;
}

Matrix& Matrix::operator=(const Matrix *rhs)
{  
  sym = rhs->sym;
  nrow = rhs->nrow;
  ncol = rhs->ncol;
  double **MC = new_dup_mat(nrow, ncol, rhs->M); // to be safe for self assignment
  if(M) delete_mat(M); 
  M = new_mat(nrow, ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = MC[j][i];
  delete_mat(MC);
  return *this;
}

Matrix& Matrix::operator+=(const Matrix &rhs)
{  
  if(sym==true && rhs.sym==false) sym= false;
  if(nrow != rhs.nrow || ncol != rhs.ncol) exit(1); //error("Trying to add matrices with different dimensions.");
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] += rhs.M[j][i];
  //printf("Matrix Compound Assigned\n");
  return *this;
}


Matrix& Matrix::operator+=(const double &x)
{  
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] += x;
  return *this;
}

Matrix& Matrix::operator*=(const double &x)
{  
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] *= x;
  return *this;
}

Matrix& Matrix::operator-(void)
{  
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = -M[j][i];
  return *this;
}


double* Matrix::operator[](const int col)
{
  return M[col];
}

const double* Matrix::operator[](const int col) const
{
  return M[col];
}

int Matrix::Rows(void) const { return nrow; }

int Matrix::Cols(void) const { return ncol; }

void Matrix::ID(void){
  sym=true;
  for(int j=0; j<ncol; j++){
    for(int i=0; i<nrow; i++){
      if(i==j) M[j][i] = 1.0;
      else M[j][i] = 0.0;
    }
  }
}

void Matrix::Zero(void){
  sym=false;  // of course, not technically true.  But safe.
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = 0.0;
}

bool Matrix::IsSym(void) const { return sym; }

void Matrix::SetSym(bool sym_in){ sym=sym_in; }

void Matrix::Scale(double val)
{
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = M[j][i]*val;
}

   
void Matrix::Shift(double val)
{
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) M[j][i] = M[j][i] + val;
}

int Matrix::Min(void) const { 
  double min = M[0][0];
  int mini = 0;
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) if(M[j][i] < min){ min=M[j][i]; mini = j*nrow+i; }
  return mini;
}

int Matrix::Max(void) const { 
  double max = M[0][0];
  int maxi = 0;
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) if(M[j][i] > max){ max=M[j][i]; maxi = j*nrow+i; }
  return maxi;
}

double Matrix::Sum(void) const {
  double sum = 0.0;
  for(int j=0; j<ncol; j++) for(int i=0; i<nrow; i++) sum += M[j][i];
  return sum;
}

void Matrix::Normalize(void) {
  double sum = Sum();
  if(sum == 0) return;
  *this *= 1.0/sum;
}


double Matrix::DetSym(void) const{
  if(!sym || nrow != ncol) error("Call of DetSym on non-sym or non-square matrix.");
  Matrix chol = Cholesky(this);

  double det = chol[0][0];
  for(int i=1; i<nrow; i++) det = det*chol[i][i];
  det = det*det;
  return det;
}


double Matrix::Trace(void) const{
  if(nrow != ncol) error("Call of Trace on non-square matrix.");
  double trace = 0.0;
  for(int i=0; i<nrow; i++) trace += M[i][i];
  return trace;
}

void Matrix::isample(int *row, int *col, void *state) const{
  int idx = -1;
  int n = ncol*nrow;
  double *probs = new_dvec(n);
  for(int i=0; i < nrow; i++) for(int j=0; j < ncol; j++)
				probs[j*nrow+i] = M[j][i];
  normalize(probs, n);
  indexsample(&idx, 1, n, probs, state);
  *col = idx/nrow;
  *row = idx%nrow;
  assert(*col < ncol);
  assert(*row < nrow);
  free(probs);
}


void Matrix::rN(const Matrix *mu, const Matrix *sig, void *state)
{      
  Zero();
  if(sig->nrow!=nrow || sig->nrow!=sig->ncol || !sig->sym || mu->ncol != 1) error("Bad rN params.");
  double *z = new_dvec(nrow);
  double **chol = new_dup_mat(nrow, nrow, sig->M);
  if(nrow==1) chol[0][0] = sqrt(chol[0][0]);
  else la_dpotrf(nrow, chol);
    
  for(int j=0;j<ncol;j++){
    for(int k=0; k<nrow; k++) z[k] = rnor(state);
    for(int i=0; i<nrow; i++){ 
      for(int k=0; k <= i; k++) M[j][i] += chol[k][i]*z[k];
      M[j][i] += mu->M[0][i];
    }
  }
  delete_mat(chol);
  free(z);
}


void Matrix::rST(double df, const Matrix *mu, const Matrix *sig, void *state){  
  if(mu->ncol != 1) error("Bad rST params.");

  Matrix none = Matrix(nrow, 1);
  rN(&none, sig, state);
  
  double v;
  
  for(int j=0;j<ncol;j++){
    v = sqrt(df/rgam(0.5*df, 0.5, state));
    for(int i=0; i<nrow; i++){
      M[j][i] = M[j][i]*v + mu->M[0][i];
    }
  }

}

void Matrix::rWSH(int nu, const Matrix *S, void *state)
{   // E = nu*S^-1
  if(nu < nrow) error("Too small df in Wishart draw.");
  Matrix A = Matrix(nrow, nu);

  Matrix none = Matrix(nrow, 1);
  Matrix Si = Inverse(S);
  A.rN(none, &Si, state);
  
  *this = A*Transpose(&A);
  sym = true;
}

void Matrix::rN(const Matrix &mu, const Matrix &sig, void *state)
{      
  Zero();
  if(sig.nrow!=nrow || sig.nrow!=sig.ncol || !sig.sym || mu.ncol != 1) error("Bad rN params.");
  double *z = new_dvec(nrow);
  double **chol = new_dup_mat(nrow, nrow, sig.M);
  if(nrow==1) chol[0][0] = sqrt(chol[0][0]);
  else la_dpotrf(nrow, chol);
    
  for(int j=0;j<ncol;j++){
    for(int k=0; k<nrow; k++) z[k] = rnor(state);
    for(int i=0; i<nrow; i++){ 
      for(int k=0; k <= i; k++) M[j][i] += chol[k][i]*z[k];
      M[j][i] += mu[0][i];
    }
  }
  delete_mat(chol);
  free(z);
}


void Matrix::rST(double df, const Matrix &mu, const Matrix &sig, void *state){  
  if(mu.ncol != 1) error("Bad rST params.");

  Matrix none = Matrix(nrow, 1);
  rN(none, sig, state);
  
  double v;
  
  for(int j=0;j<ncol;j++){
    v = sqrt(df/rgam(0.5*df, 0.5, state));
    for(int i=0; i<nrow; i++){
      M[j][i] = M[j][i]*v + mu[0][i];
    }
  }

}

void Matrix::rWSH(int nu, const Matrix &S, void *state)
{   // E = nu*S^-1
  if(nu < nrow) error("Too small df in Wishart draw.");
  Matrix A = Matrix(nrow, nu);

  Matrix none = Matrix(nrow, 1);
  Matrix Si = Inverse(S);
  A.rN(none, Si, state);
  
  *this = A*Transpose(A);
  sym = true;
}

double Matrix::ldN(const Matrix *mu, const Matrix *sig) const
{   
  
  if(mu->ncol != 1) error("Bad ldN params.");
 
  /*  get the normalizing constant */
  double det = sig->DetSym();
  double lconst = -0.5*( ((double) nrow)*log(2.0*PI) + log(det) );

  /* calculate the lhd */

  double loglhd = lconst*((double) ncol);

  for(int j=0; j<ncol; j++){
    Matrix diff = Matrix(nrow, 1, M[j]) - *mu;
    Matrix R = Transpose(&diff)*Inverse(sig)*diff;
    loglhd += -0.5*R[0][0];
  }

  return loglhd;
}

double Matrix::ldN(const Matrix &mu, const Matrix &sig) const
{   
  
  if(mu.ncol != 1) error("Bad ldN params.");
 
  /*  get the normalizing constant */
  double det = sig.DetSym();
  double lconst = -0.5*( ((double) nrow)*log(2.0*PI) + log(det) );

  /* calculate the lhd */

  double loglhd = lconst*((double) ncol);

  for(int j=0; j<ncol; j++){
    Matrix diff = Matrix(nrow, 1, M[j]) - mu;
    Matrix R = Transpose(diff)*Inverse(sig)*diff;
    loglhd += -0.5*R[0][0];
  }

  return loglhd;
}


double Matrix::ldST(double df, const Matrix *mu, const Matrix *sig) const
{  
  if(mu->Cols() != 1) error("Bad ldST params.");
  
  if(nrow==1){
    double ld = 0.0;
    for(int j=0; j<ncol; j++) 
      ld += log(dt( (M[j][0] - mu->M[0][0])/sqrt(sig->M[0][0]), df, 0)/sqrt(sig->M[0][0]));
    return ld; }

  double dim = (double) nrow;

  /*  get the normalizing constant */
    
  double lconst = -0.5*(log(sig->DetSym())+dim*log(df*PI));
  lconst += lgammafn( (df + dim)*0.5 ) - lgammafn( 0.5*df );

  /*  the rest */

  double loglhd = lconst*((double) ncol);

  for(int j=0; j<ncol; j++){
    Matrix diff = Matrix(nrow, 1, M[j]) - *mu;
    Matrix R = Transpose(&diff)*Inverse(sig)*diff;
    loglhd += -0.5*(df+dim)*log( 1.0 + R[0][0]/df );
  }

  return loglhd;
}


double Matrix::ldST(double df, const Matrix &mu, const Matrix &sig) const
{  
  if(mu.ncol != 1) error("Bad ldST params.");
  
  if(nrow==1){
    double ld = 0.0;
    for(int j=0; j<ncol; j++) 
      ld += log(dt( (M[j][0] - mu.M[0][0])/sqrt(sig.M[0][0]), df, 0)/sqrt(sig.M[0][0]));
    return ld; }

  double dim = (double) nrow;

  /*  get the normalizing constant */
    
  double lconst = -0.5*(log(sig.DetSym())+dim*log(df*PI));
  lconst += lgammafn( (df + dim)*0.5 ) - lgammafn( 0.5*df );

  /*  the rest */

  double loglhd = lconst*((double) ncol);

  for(int j=0; j<ncol; j++){
    Matrix diff = Matrix(nrow, 1, M[j]) - mu;
    Matrix R = Transpose(diff)*Inverse(sig)*diff;
    loglhd += -0.5*(df+dim)*log( 1.0 + R[0][0]/df );
  }

  return loglhd;
}


double Matrix::ldWSH(int nu, const Matrix *B) const //E[W] = nu*B^-1
{
  double d = (double) nrow;
  double v = (double) nu;

  double lconst = v*log(B->DetSym())-d*(d-1)*0.25*log(PI);
  for(double i=1.0; i<=d; i++) lconst += -lgammafn( (2.0*v +1.0 - i)*0.5 );

  double loglhd = (2.0*v - d - 1)*0.5*log(this->DetSym());
  loglhd += lconst;

  Matrix BW = (*B)*(*this);
  loglhd += -BW.Trace();
  
  return loglhd;
}

double Matrix::ldWSH(int nu, const Matrix &B) const //E[W] = nu*B^-1
{
  double d = (double) nrow;
  double v = (double) nu;

  double lconst = v*log(B.DetSym())-d*(d-1)*0.25*log(PI);
  for(double i=1.0; i<=d; i++) lconst += -lgammafn( (2.0*v +1.0 - i)*0.5 );

  double loglhd = (2.0*v - d - 1)*0.5*log(this->DetSym());
  loglhd += lconst;

  Matrix BW = B*(*this);
  loglhd += -BW.Trace();
  
  return loglhd;
}


/*
 * return an integer of length (*len) wih indexes into V[][col] which
 * satisfy the relation "V op val" where op is one of LT(<) GT(>)
 * EQ(==) LEQ(<=) GEQ(>=) NE(!=)
 */

int* Matrix::find_col(int var, FIND_OP op, double val, int *len){
  
  int i,j;
  int *tf = new_ivec(ncol);
  
  assert(var < nrow);

  (*len) = 0;
  switch (op) {
  case GT:  
    for(i=0; i<ncol; i++) {
      if(M[i][var] >  val) tf[i] = 1; 
      else tf[i] = 0; 
      if(tf[i] == 1) (*len)++;
    }
    break;
  case GEQ: 
    for(i=0; i<ncol; i++) {
      if(M[i][var] >= val) tf[i] = 1; 
      else tf[i] = 0; 
      if(tf[i] == 1) (*len)++;
    }
    break;
  case EQ:  
    for(i=0; i<ncol; i++) {
      if(M[i][var] == val) tf[i] = 1; 
      else tf[i] = 0; 
      if(tf[i] == 1) (*len)++;
    }
    break;
  case LEQ: 
    for(i=0; i<ncol; i++) {
      if(M[i][var] <= val) tf[i] = 1; 
      else tf[i] = 0; 
      if(tf[i] == 1) (*len)++;
    }
    break;
  case LT:  
    for(i=0; i<ncol; i++) {
      if(M[i][var] <  val) tf[i] = 1; 
      else tf[i] = 0; 
      if(tf[i] == 1) (*len)++;
    }
    break;
  case NE:  
    for(i=0; i<ncol; i++) {
      if(M[i][var] != val) tf[i] = 1; 
      else tf[i] = 0; 
      if(tf[i] == 1) (*len)++;
    }
    break;
  default: error("OP not supported");
  }
  
  int *found;
  if(*len == 0) found = NULL;
  else {
    found = new_ivec(*len);
    for(i=0,j=0; i<ncol; i++) {
      if(tf[i]) {
	found[j] = i;
	j++;
      }
    }
  }
  
  free(tf);
  return found;

}



double* Matrix::Mean(int index){

  Matrix tmp;
  if(index==0)  tmp = Transpose(*this);
  else tmp = Matrix(*this);
  
  double* ret = new_dzero(tmp.Cols());

  for(int j=0; j<tmp.Cols(); j++){
    for(int i=0; i<tmp.Rows(); i++) ret[j] += tmp[j][i];
    ret[j] = ret[j]/((double) tmp.Rows());
  }
  return ret;
}

/******************** friends ************************/


Matrix cbind(const Matrix& lhs, const Matrix& rhs){

  if(lhs.Rows()!=rhs.Rows()) error("Bad dimensions in cbind.");
  int nrow = lhs.Rows();
  int ncol = lhs.Cols() + rhs.Cols();

  
  double *vals = new_dvec(nrow*ncol);
  for(int j=0; j<lhs.Cols(); j++) for(int i=0; i<nrow; i++) vals[j*nrow+i] = lhs[j][i];
  for(int j=0; j<rhs.Cols(); j++) for(int i=0; i<nrow; i++) vals[(j+lhs.Cols())*nrow+i] = rhs[j][i];

  Matrix ret = Matrix(nrow, ncol, vals);
  free(vals);
  return ret;

}


Matrix cbind(const Matrix *lhs, const Matrix *rhs){

  if(lhs->Rows()!=rhs->Rows()) error("Bad dimensions in cbind.");
  int nrow = lhs->Rows();
  int ncol = lhs->Cols() + rhs->Cols();

  
  double *vals = new_dvec(nrow*ncol);
  for(int j=0; j<lhs->Cols(); j++) for(int i=0; i<nrow; i++) vals[j*nrow+i] = lhs->M[j][i];
  for(int j=0; j<rhs->Cols(); j++) for(int i=0; i<nrow; i++) vals[(j+lhs->Cols())*nrow+i] = rhs->M[j][i];

  Matrix ret = Matrix(nrow, ncol, vals);
  free(vals);
  return ret;

}

Matrix rbind(const Matrix& lhs, const Matrix& rhs){

  if(lhs.Cols()!=rhs.Cols()) error("Bad dimensions in rbind.");
  int ncol = lhs.Cols();
  int nrow = lhs.Rows() + rhs.Rows();

  
  double *vals = new_dvec(nrow*ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<lhs.Rows(); i++) vals[j*nrow+i] = lhs[j][i];
  for(int j=0; j<ncol; j++) for(int i=0; i<rhs.Rows(); i++)  vals[j*nrow+lhs.Rows() + i] = rhs[j][i];

  Matrix ret = Matrix(nrow, ncol, vals);
  free(vals);
  return ret;

}


Matrix rbind(const Matrix *lhs, const Matrix *rhs){

  if(lhs->Cols()!=rhs->Cols()) error("Bad dimensions in rbind.");
  int ncol = lhs->Cols();
  int nrow = lhs->Rows() + rhs->Rows();

  
  double *vals = new_dvec(nrow*ncol);
  for(int j=0; j<ncol; j++) for(int i=0; i<lhs->Rows(); i++) vals[j*nrow+i] = lhs->M[j][i];
  for(int j=0; j<ncol; j++) for(int i=0; i<rhs->Rows(); i++)  vals[j*nrow+lhs->Rows() + i] = rhs->M[j][i];

  Matrix ret = Matrix(nrow, ncol, vals);
  free(vals);
  return ret;

}


Matrix operator*(const Matrix& lhs, const Matrix& rhs){   
  if(lhs.Cols()!=rhs.Rows()) error("Bad dimensions in Matrix multiply.");

  Matrix ret = Matrix(lhs.Rows(), rhs.Cols());

  Matrix L = Matrix(lhs);
  Matrix R = Matrix(rhs);

  if(lhs.Rows() == 1 && rhs.Rows() == 1 && lhs.Cols() == 1 && rhs.Cols() == 1) ret.M[0][0] = (L.M[0][0])*(R.M[0][0]);
  else if(L.IsSym()) la_dsymm(true, L.Rows(), L.Cols(), 
			   R.Rows(), R.Cols(), ret.Rows(), ret.Cols(),
			   L.M, R.M, ret.M, 1.0, 0.0);
  else if(R.IsSym()) la_dsymm(false, R.Rows(), R.Cols(), 
			   L.Rows(), L.Cols(), ret.Rows(), ret.Cols(),
			   R.M, L.M, ret.M, 1.0, 0.0);
  else la_dgemm(false, false, L.Rows(), L.Cols(), 
			   R.Rows(), R.Cols(), ret.Rows(), ret.Cols(),
			   L.M, R.M, ret.M, 1.0, 0.0);
  
  return ret;
}


Matrix Inverse(const Matrix *mat){

  if(mat->Rows()!=mat->Cols()) error("Attempt to invert non-square matrix.");

  Matrix mati = Matrix(mat);
  if(mat->Rows() == 1){ mati.M[0][0] = 1.0/mat->M[0][0]; mati.SetSym(true); return mati; }

  double **I = new_id_mat(mat->Rows());

  if(mat->IsSym()) la_dposv(mat->Rows(), mat->Cols(), mati.M, I);
  else la_dgesv(mat->Rows(), mat->Cols(), mati.M, I);

  copy_mat(mat->Rows(), mat->Cols(), mati.M, I);
  delete_mat(I);

  return mati;
}


Matrix Inverse(const Matrix& mat){

  if(mat.Rows()!=mat.Cols()) error("Attempt to invert non-square matrix.");

  Matrix mati = Matrix(mat);
  if(mat.Rows() == 1){ mati.M[0][0] = 1.0/mat.M[0][0]; mati.SetSym(true); return mati; }

  double **I = new_id_mat(mat.Rows());

  if(mat.IsSym()) la_dposv(mat.Rows(), mat.Cols(), mati.M, I);
  else la_dgesv(mat.Rows(), mat.Cols(), mati.M, I);

  copy_mat(mat.Rows(), mat.Cols(), mati.M, I);
  delete_mat(I);

  return mati;
}

// returns the lower-tri cholesky ;  zeros out the strict upper triangle.

Matrix Cholesky(const Matrix& mat){
  if(!mat.IsSym() || mat.Rows()!=mat.Cols()) error("Cholesky attempt on non-sym or non-square matrix.");

  Matrix chol = Matrix(mat);
  if(mat.Rows() == 1){ chol.M[0][0] = sqrt(mat.M[0][0]); return mat; }

  la_dpotrf(chol.Rows(), chol.M);
  for(int i=0; i<chol.Rows(); i++) for(int j=i+1; j<chol.Cols(); j++) chol[j][i] = 0.0;
  chol.SetSym(false);

  return chol;
}


// returns the lower-tri cholesky ;  zeros out the strict upper triangle.

Matrix Cholesky(const Matrix *mat){
  if(!mat->IsSym() || mat->Rows()!=mat->Cols()) error("Cholesky attempt on non-sym or non-square matrix.");

  Matrix chol = Matrix(mat);
  if(mat->Rows() == 1){ chol.M[0][0] = sqrt(mat->M[0][0]); return mat; }

  la_dpotrf(chol.Rows(), chol.M);
  for(int i=0; i<chol.Rows(); i++) for(int j=i+1; j<chol.Cols(); j++) chol[j][i] = 0.0;
  chol.SetSym(false);

  return chol;
}


/******************** non-member functions **********************/	 


Matrix Logit(const Matrix& mat){
  Matrix lmat = Matrix(mat.Rows(), mat.Cols());
  for(int i=0; i<mat.Rows(); i++) for(int j=0; j<mat.Cols(); j++) lmat[j][i] = logit(mat[j][i]);
  if(mat.IsSym()) lmat.SetSym(true); 
  return lmat;
}


Matrix Expit(const Matrix& mat){
  Matrix emat = Matrix(mat.Rows(),mat.Cols());
  for(int i=0; i<mat.Rows(); i++) for(int j=0; j<mat.Cols(); j++) emat[j][i] = expit(mat[j][i]);
  if(mat.IsSym()) emat.SetSym(true); 
  return emat;
}

Matrix Transpose(const Matrix& mat){
  Matrix Tmat = Matrix(mat.Cols(), mat.Rows());
  for(int i=0; i<mat.Rows(); i++) for(int j=0; j<mat.Cols(); j++) Tmat[i][j] += mat[j][i];
  if(mat.IsSym()) Tmat.SetSym(true); 
  return Tmat;
}


Matrix Logit(const Matrix *mat){
  Matrix lmat = Matrix(mat->Rows(), mat->Cols());
  for(int i=0; i<mat->Rows(); i++) for(int j=0; j<mat->Cols(); j++) lmat[j][i] = logit(mat->M[j][i]);
  if(mat->IsSym()) lmat.SetSym(true); 
  return lmat;
}


Matrix Expit(const Matrix *mat){
  Matrix emat = Matrix(mat->Rows(),mat->Cols());
  for(int i=0; i<mat->Rows(); i++) for(int j=0; j<mat->Cols(); j++) emat[j][i] = expit(mat->M[j][i]);
  if(mat->IsSym()) emat.SetSym(true); 
  return emat;
}

Matrix Transpose(const Matrix *mat){
  Matrix Tmat = Matrix(mat->Cols(), mat->Rows());
  for(int i=0; i<mat->Rows(); i++) for(int j=0; j<mat->Cols(); j++) Tmat[i][j] += mat->M[j][i];
  if(mat->IsSym()) Tmat.SetSym(true); 
  return Tmat;
}

ostream& operator<<(ostream& out, const Matrix& mat){
  if(mat.Rows()==0 || mat.Cols()==0){ out << "NULL" << endl; return out; }
  out << endl;
  for(int i=0; i<mat.Rows(); i++){
    for(int j=0; j<mat.Cols(); j++) out << mat[j][i] << " ";
    out << endl;
  }
  return out;
}   

Matrix operator+(const Matrix& lhs, const Matrix& rhs){  
  Matrix ret(lhs);
  ret += rhs;
  return ret;
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs){  
  Matrix ret(rhs);
  -ret += lhs;
  return ret;
}


Matrix operator+(const Matrix& lhs, const double& x){
  Matrix ret(lhs);
  ret += x;
  return ret;
}

Matrix operator+(const double& x, const Matrix& rhs){ return rhs+x; }

Matrix operator-(const Matrix& lhs, const double& x){
  Matrix ret(lhs);
  ret += -x;
  return ret;
}

Matrix operator-(const double& x, const Matrix& rhs){ 
   Matrix ret(rhs);
  -ret += x;
  return ret;
}

Matrix operator*(const Matrix& lhs, const double& x){
  Matrix ret(lhs);
  ret *= x;
  return ret;
}

Matrix operator*(const double& x, const Matrix& rhs){ return rhs*x; }

Matrix operator/(const Matrix& lhs, const double& x){
  Matrix ret(lhs);
  ret *= (1.0/x);
  return ret;
}
