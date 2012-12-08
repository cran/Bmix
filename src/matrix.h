/* A basic data class.  Written by Matt Taddy */

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>


typedef enum FIND_OP {LT=101, LEQ=102, EQ=103, GEQ=104, GT=105, NE=106} FIND_OP;

using namespace std;


class Matrix
{
 private:

  bool sym;
  int nrow;
  int ncol;

  /* the actual values are stored (in M) column major, 
     even though the interface is row major */

  double **M; 

 
 public:
  /* The matrix is filled column major, as in R */
  Matrix();
  Matrix(int nrow_in, int ncol_in);
  Matrix(int nrow_in, int ncol_in, double val);
  Matrix(int nrow_in, int ncol_in, const double *vals);
  Matrix(int nrow_in, int ncol_in, const double *vals, bool sym_in); 
  Matrix(int nrow_in, int ncol_in, const Matrix &mat, int *rows, int *cols);
  Matrix(int nrow_in, int ncol_in, const Matrix *mat, int *rows, int *cols);
  Matrix(const Matrix& mat);
  Matrix(const Matrix *mat);
  ~Matrix();
  

  Matrix& operator=(const Matrix& rhs);
  Matrix& operator=(const Matrix *rhs);
  Matrix& operator+=(const Matrix& rhs);
  Matrix& operator+=(const Matrix *rhs);
  Matrix& operator+=(const double& x);
  Matrix& operator*=(const double& x);
  Matrix& operator-(void);
  double* operator[](const int index);
  const double* operator[](const int index) const;

  int Rows() const;
  int Cols() const;

  void Zero();
  void ID();

  bool IsSym() const;
  void SetSym(bool sym_in);

  void Scale(double val); // multiply M by val
  void Shift(double val); // add val to M
  
  int Min() const;
  int Max() const;
  double* Mean(int index);  
  double Sum() const;
  void Normalize();
  double DetSym() const;
  double Trace() const;
  void isample(int *row, int *col, void *state) const;

  void rN(const Matrix *mu, const Matrix *sig, void *state);
  void rST(double df, const Matrix *mu, const Matrix *sig, void *state);
  void rWSH(int nu, const Matrix *S, void *state);

  void rN(const Matrix &mu, const Matrix &sig, void *state);
  void rST(double df, const Matrix &mu, const Matrix &sig, void *state);
  void rWSH(int nu, const Matrix &S, void *state);

  double ldN(const Matrix *mu, const Matrix *sig) const;
  double ldST(double df, const Matrix *mu, const Matrix *sig) const;
  double ldWSH(int nu, const Matrix *S) const;

  double ldN(const Matrix &mu, const Matrix &sig) const;
  double ldST(double df, const Matrix &mu, const Matrix &sig) const;
  double ldWSH(int nu, const Matrix &S) const;

  int* find_col(int var, FIND_OP op, double val, int *nc);
	    
  friend Matrix Inverse(const Matrix& mat);  
  friend Matrix Cholesky(const Matrix& mat); // lower-tri 
  friend Matrix Inverse(const Matrix *mat);  
  friend Matrix Cholesky(const Matrix *mat); // lower-tri 

  friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
  friend Matrix cbind(const Matrix& lhs, const Matrix& rhs);
  friend Matrix cbind(const Matrix *lhs, const Matrix *rhs);
  friend Matrix rbind(const Matrix& lhs, const Matrix& rhs);
  friend Matrix rbind(const Matrix *lhs, const Matrix *rhs);

  friend Matrix Logit(const Matrix *mat);
  friend Matrix Expit(const Matrix *mat);
  friend Matrix Transpose(const Matrix *mat);
  friend Matrix Logit(const Matrix& mat);
  friend Matrix Expit(const Matrix& mat);
  friend Matrix Transpose(const Matrix& mat);

};

ostream& operator<<(ostream& out, const Matrix& mat);
Matrix operator+(const Matrix& lhs, const Matrix& rhs);
Matrix operator-(const Matrix& lhs, const Matrix& rhs);
Matrix operator+(const Matrix& lhs, const double& x);
Matrix operator+(const double& x, const Matrix& rhs);
Matrix operator-(const Matrix& lhs, const double& x);
Matrix operator-(const double& x, const Matrix& rhs);
Matrix operator*(const Matrix& lhs, const double& x);
Matrix operator*(const double& x, const Matrix& rhs);
Matrix operator/(const Matrix& lhs, const double& x);



#endif
  

 
