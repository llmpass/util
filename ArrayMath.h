#ifndef _ARRAYMATH_H
#define _ARRAYMATH_H

namespace util {

  #define PI 3.1415926535
  #define maxf  1e30
  #define minf -1e30

  ////////////////////////////////////////////////////////////////////////////
  // min max
  float minv(float* x, int n, int& ind); 
  float maxv(float* x, int n, int& ind);
  float minv(float** x, int n1, int n2, int& ind1, int& ind2); 
  float maxv(float** x, int n1, int n2, int& ind1, int& ind2);
  double minv(double* x, int n, int& ind); 
  double maxv(double* x, int n, int& ind);
  double minv(double** x, int n1, int n2, int& ind1, int& ind2); 
  double maxv(double** x, int n1, int n2, int& ind1, int& ind2);

  ////////////////////////////////////////////////////////////////////////////
  // copy 
  void copy(float* x, float* y, const int n); 
  float* copy(float* x, const int n);
  void copy(float** x, float** y, const int n1, const int n2);
  float** copy(float **x, const int n1, const int n2);
  void copy(float*** x, float*** y, const int n1, const int n2, const int n3);
  float*** copy(float ***x, const int n1, const int n2, const int n3);
  void copy(double* x, double* y, const int n);
  double* copy(double* x, const int n);
  void copy(double** x, double** y, const int n1, const int n2);
  double** copy(double **x, const int n1, const int n2);
  void copy(double*** x, double*** y, const int n1, const int n2, 
    const int n3);
  double*** copy(double ***x, const int n1, const int n2, const int n3);
  void copy(const int nn, float* x, float* y);
  void copy(const int nn1, const int nn2, float** x, float** y);
  void copy(const int nn1, const int nn2, const int nn3, 
    float*** x, float*** y);
  void ccopy(int n1, float* cx, float* cy);
  void ccopy(int n1, int n2, float** cx, float** cy);
  void ccopy(int n1, int n2, int n3, float*** cx, float*** cy);
  void fill(float a, float* x, int n);
  void fill(float a, float** x, int n1, int n2);
  float* fillfloat(float ra, int n1);
  void fill(float a, float*** x, int n1, int n2, int n3);
  float** fillfloat(float ra, int n1, int n2);
  float*** fillfloat(float ra, int n1, int n2, int n3);
  void fill(int a, int* x, int n);
  void fill(int a, int** x, int n1, int n2);
  void fill(int a, int*** x, int n1, int n2, int n3);
  void fill(short a, short* x, int n);
  void fill(short a, short** x, int n1, int n2);
  void fill(short a, short*** x, int n1, int n2, int n3);
  void fill(double a, double* x, int n);
  void fill(double a, double **x, int n1, int n2);
  void fill(double a, double ***x, int n1, int n2, int n3);
  void zero(int* x, int n);
  void zero(int** x, int n1, int n2);
  void zero(int*** x, int n1, int n2, int n3);
  void zero(float* x, int n);
  void zero(float** x, int n1, int n2);
  void zero(float*** x, int n1, int n2, int n3);
  void zero(double* x, int n);
  void zero(double** x, int n1, int n2);
  void zero(double*** x, int n1, int n2, int n3);
  int* zeroint(int n);
  int** zeroint(int n1, int n2);
  int*** zeroint(int n1, int n2, int n3);
  float* zerofloat(int n);
  float** zerofloat(int n1, int n2);
  float*** zerofloat(int n1, int n2, int n3);
  double* zerodouble(int n);
  double** zerodouble(int n1, int n2);
  double*** zerodouble(int n1, int n2, int n3);
  int binarySearch(int* a, int alength, int x, int i);
  int binarySearch(int* a, int alength, int x);
  void div(float ra, float *ry, float *rz, int n);
  void div(float ra, float **ry, float **rz, int n1, int n2);
  void div(float ra, float ***ry, float ***rz, int n1, int n2, int n3);
  void mul(float ra, float *ry, float *rz, int n);
  void mul(float ra, float **ry, float **rz, int n1, int n2);
  float* mul(float ra, float *rx, int n);
  void add(float ra, float *ry, float *rz, int n);
  void add(float ra, float **ry, float **rz, int n1, int n2);
  void add(float ra, float ***ry, float ***rz, int n1, int n2, int n3);
  void add(float* rx, float *ry, float *rz, int n);
  void add(float** rx, float **ry, float **rz, int n1, int n2);
  void add(float*** rx, float ***ry, float ***rz, int n1, int n2, int n3);
  float* add(float* x, float* y, int n);
  float** add(float** x, float** y, int n1, int n2);
  float*** add(float*** x, float*** y, int n1, int n2, int n3);
  void sub(float ra, float *ry, float *rz, int n);
  void sub(float ra, float **ry, float **rz, int n1, int n2);
  void sub(float ra, float ***ry, float ***rz, int n1, int n2, int n3);
  float* sub(float ra, float* ry, int n);
  float** sub(float ra, float** ry, int n1, int n2);
  float*** sub(float ra, float*** ry, int n1, int n2, int n3);
  void sub(float* rx, float *ry, float *rz, int n);
  void sub(float** rx, float **ry, float **rz, int n1, int n2);
  void sub(float*** rx, float ***ry, float ***rz, int n1, int n2, int n3);
  float* sub(float* x, float* y, int n);
  float** sub(float** x, float** y, int n1, int n2);
  float*** sub(float*** x, float*** y, int n1, int n2, int n3);
  void sub(double ra, double *ry, double *rz, int n);
  void sub(double ra, double **ry, double **rz, int n1, int n2);
  void sub(double ra, double ***ry, double ***rz, int n1, int n2, int n3);
  double* sub(double ra, double* ry, int n);
  double** sub(double ra, double** ry, int n1, int n2);
  double*** sub(double ra, double*** ry, int n1, int n2, int n3);
  void sub(double* rx, double *ry, double *rz, int n);
  void sub(double** rx, double **ry, double **rz, int n1, int n2);
  void sub(double*** rx, double ***ry, double ***rz, int n1, int n2, int n3);
  double* sub(double* x, double* y, int n);
  double** sub(double** x, double** y, int n1, int n2);
  double*** sub(double*** x, double*** y, int n1, int n2, int n3);
  double* float2double(float* f, int n);
  double** float2double(float** f, int n1, int n2);
  double*** float2double(float*** f, int n1, int n2, int n3);
  float* double2float(double* f, int n);
  float** double2float(double** f, int n1, int n2);
  float*** double2float(double*** f, int n1, int n2, int n3);
  void mulmat(float** a, float** b, float** c, int n, int m, int l);
  void mulmat(float** a, float *v, float *u, int n, int m);
  void mulmat(double **a, double **b, double **c, int n, int m, int l);
  void mulmat(double **a, double *v, double *u, int n, int m);
  float sum(float* x, int n);
  float sum(float** x, int n1, int n2);
  float sum(float*** x, int n1, int n2, int n3);
  double sum(double* x, int n);
  double sum(double** x, int n1, int n2);
  double sum(double*** x, int n1, int n2, int n3);
  void float2double(float* x, double* y, int n1);
  void float2double(float** x, double** y, int n1, int n2); 
  void float2double(float*** x, double*** y, int n1, int n2, int n3);
  void double2float(double* x, float* y, int n1);
  void double2float(double** x, float** y, int n1, int n2);
  void double2float(double*** x, float*** y, int n1, int n2, int n3);
  double correlation(double** f, double** g, int n1, int n2);
}

#endif
