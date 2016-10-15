#include "ArrayMath.h"
#include <cmath>
#include <iostream>
using namespace std;

float util::minv(float* x, int n, int& ind) {
  float m = maxf;
  for (int i=0; i<n; ++i) 
    if (m>x[i]) {m = x[i]; ind = i;}
  return m;
}
float util::maxv(float* x, int n, int& ind) {
  float m = minf;
  for (int i=0; i<n; ++i) 
    if (m<x[i]) {m = x[i]; ind = i;}
  return m;
}
float util::minv(float** x, int n1, int n2, int& ind1, int& ind2) {
  float m = maxf;
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) 
      if (m>x[i2][i1]) {m = x[i2][i1]; ind1 = i1; ind2 = i2;}
  return m;
}
float util::maxv(float** x, int n1, int n2, int& ind1, int& ind2) {
  float m = minf;
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) 
      if (m<x[i2][i1]) {m = x[i2][i1]; ind1 = i1; ind2 = i2;}
  return m;
}
double util::minv(double* x, int n, int& ind) {
  double m = maxf;
  for (int i=0; i<n; ++i) 
    if (m>x[i]) {m = x[i]; ind = i;}
  return m;
}
double util::maxv(double* x, int n, int& ind) {
  double m = minf;
  for (int i=0; i<n; ++i) 
    if (m<x[i]) {m = x[i]; ind = i;}
  return m;
}
double util::minv(double** x, int n1, int n2, int& ind1, int& ind2) {
  double m = maxf;
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) 
      if (m>x[i2][i1]) {m = x[i2][i1]; ind1 = i1; ind2 = i2;}
  return m;
}
double util::maxv(double** x, int n1, int n2, int& ind1, int& ind2) {
  double m = minf;
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) 
      if (m<x[i2][i1]) {m = x[i2][i1]; ind1 = i1; ind2 = i2;}
  return m;
}

void util::copy(float* x, float* y, const int n) {
  for (int i=0; i<n; ++i) y[i]=x[i];
}
float* util::copy(float* x, const int n) {
  float* y = new float[n];
  copy(x,y,n);
  return y;
}
void util::copy(float** x, float** y, const int n1, const int n2) {
  for (int i2=0; i2<n2; ++i2) copy(x[i2],y[i2],n1);
}
float** util::copy(float **x, const int n1, const int n2) {
  float** y = new float*[n2];
  for (int i2=0; i2<n2; ++i2) y[i2] = new float[n1];
  copy(x,y,n1,n2);
  return y;
}
void util::copy(float*** x, float*** y, 
  const int n1, const int n2, const int n3) {
  for (int i3=0; i3<n3; ++i3) copy(x[i3],y[i3],n1,n2);
}
float*** util::copy(float ***x, const int n1, const int n2, const int n3) {
  float*** y = new float**[n3];
  for (int i3=0; i3<n3; ++i3) {
    y[i3] = new float*[n2];
    for (int i2=0; i2<n2; ++i2) y[i3][i2] = new float[n1];
  }
  copy(x,y,n1,n2,n3);
  return y;
}
void util::copy(double* x, double* y, const int n) {
  for (int i=0; i<n; ++i) y[i]=x[i];
}
double* util::copy(double* x, const int n) {
  double* y = new double[n];
  copy(x,y,n);
  return y;
}
void util::copy(double** x, double** y, const int n1, const int n2) {
  for (int i2=0; i2<n2; ++i2) copy(x[i2],y[i2],n1);
}
double** util::copy(double **x, const int n1, const int n2) {
  double** y = new double*[n2];
  for (int i2=0; i2<n2; ++i2) y[i2] = new double[n1];
  copy(x,y,n1,n2);
  return y;
}
void util::copy(double*** x, double*** y, 
  const int n1, const int n2, const int n3) {
    for (int i3=0; i3<n3; ++i3) copy(x[i3],y[i3],n1,n2);
  }
double*** util::copy(double ***x, const int n1, const int n2, const int n3) {
  double*** y = new double**[n3];
  for (int i3=0; i3<n3; ++i3) {
    y[i3] = new double*[n2];
    for (int i2=0; i2<n2; ++i2) y[i3][i2] = new double[n1];
  }
  copy(x,y,n1,n2,n3);
  return y;
}
 
/**
 * Copies elements from one specified array to another.
 * @param nn number of elements to copy in 1st dimension.
 * @param x source array.
 * @param y destination array.
 */
void util::copy(const int nn, float* x, float* y) {
  for (int i1=0; i1<nn; ++i1) y[i1] = x[i1];
}
void util::copy(const int nn1, const int nn2, float** x, float** y) {
  for (int i2=0; i2<nn2; ++i2) 
    for (int i1=0; i1<nn1; ++i1) y[i2][i1] = x[i2][i1];
}
void util::copy(const int nn1, const int nn2, const int nn3, 
  float*** x, float*** y) {
  for (int i3=0; i3<nn3; ++i3) 
    for (int i2=0; i2<nn2; ++i2) 
      for (int i1=0; i1<nn1; ++i1) y[i3][i2][i1] = x[i3][i2][i1];
}
void util::ccopy(int n1, float* cx, float* cy) {
  copy(2*n1,cx,cy);
}
void util::ccopy(int n1, int n2, float** cx, float** cy) {
  for (int i2=0; i2<n2; ++i2)
    ccopy(n1,cx[i2],cy[i2]);
}
void util::ccopy(int n1, int n2, int n3, float*** cx, float*** cy) {
  for (int i3=0; i3<n3; ++i3)
    ccopy(n1,n2,cx[i3],cy[i3]);
}

  
////////////////////////////////////////////////////////////////////////////
// fill
  
void util::fill(float a, float* x, int n) {
  for (int i=0; i<n; ++i) x[i]=a;
}
void util::fill(float a, float** x, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) x[i2][i1]=a;
}
float* util::fillfloat(float ra, int n1) {
  float* rx = new float[n1];
  fill(ra,rx,n1);
  return rx;
}
void util::fill(float a, float*** x, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) x[i3][i2][i1]=a;
}

float** util::fillfloat(float ra, int n1, int n2) {
  float** rx = new float*[n2];
  for (int i2=0; i2<n2; ++i2) rx[i2] = fillfloat(ra,n1);
  return rx;
}
float*** util::fillfloat(float ra, int n1, int n2, int n3) {
  float*** rx = new float**[n3];
  for (int i3=0; i3<n3; ++i3) rx[i3] = fillfloat(ra,n1,n2);
  return rx;
}
  
void util::fill(short a, short* x, int n) {
  for (int i=0; i<n; ++i) x[i]=a;
}
void util::fill(short a, short** x, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) x[i2][i1]=a;
}
void util::fill(short a, short*** x, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) x[i3][i2][i1]=a;
}

void util::fill(int a, int* x, int n) {
  for (int i=0; i<n; ++i) x[i]=a;
}
void util::fill(int a, int** x, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) x[i2][i1]=a;
}
void util::fill(int a, int*** x, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) x[i3][i2][i1]=a;
}


void util::fill(double a, double* x, int n) {
  for (int i=0; i<n; ++i) x[i]=a;
}
void util::fill(double a, double **x, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) x[i2][i1]=a;
}
void util::fill(double a, double ***x, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) x[i3][i2][i1]=a;
}
  
////////////////////////////////////////////////////////////////////////////
// zero
  
void util::zero(int* x, int n) {fill(0,x,n);}
void util::zero(int** x, int n1, int n2) {fill(0,x,n1,n2);} 
void util::zero(int*** x, int n1, int n2, int n3) {fill(0,x,n1,n2,n3);}
void util::zero(float* x, int n) {fill(0,x,n);}
void util::zero(float** x, int n1, int n2) {fill(0,x,n1,n2);} 
void util::zero(float*** x, int n1, int n2, int n3) {fill(0,x,n1,n2,n3);}
void util::zero(double* x, int n) {fill(0,x,n);}
void util::zero(double** x, int n1, int n2) {fill(0,x,n1,n2);}
void util::zero(double*** x, int n1, int n2, int n3) {fill(0,x,n1,n2,n3);}

int* util::zeroint(int n) {
  int *x = new int[n];
  zero(x,n);
  return x;
}
int** util::zeroint(int n1, int n2) {
  int **x = new int*[n2];
  for (int i2=0; i2<n2; ++i2) x[i2]=zeroint(n1);
  return x;
}
int*** util::zeroint(int n1, int n2, int n3) {
  int ***x = new int**[n3];
  for (int i3=0; i3<n3; ++i3) x[i3]=zeroint(n1,n2);
  return x;
}

float* util::zerofloat(int n) {
  float *x = new float[n];
  zero(x,n);
  return x;
}
float** util::zerofloat(int n1, int n2) {
  float **x = new float*[n2];
  for (int i2=0; i2<n2; ++i2) x[i2]=zerofloat(n1);
  return x;
}
float*** util::zerofloat(int n1, int n2, int n3) {
  float ***x = new float**[n3];
  for (int i3=0; i3<n3; ++i3) x[i3]=zerofloat(n1,n2);
  return x;
}

double* util::zerodouble(int n) {
  double *x = new double[n];
  zero(x,n);
  return x;
}
double** util::zerodouble(int n1, int n2) {
  double **x = new double*[n2];
  for (int i2=0; i2<n2; ++i2) x[i2]=zerodouble(n1);
  return x;
}
double*** util::zerodouble(int n1, int n2, int n3) {
  double ***x = new double**[n3];
  for (int i3=0; i3<n3; ++i3) x[i3]=zerodouble(n1,n2);
  return x;
}

////////////////////////////////////////////////////////////////////////////
// binary search

/**
 * Performs a binary search in a monotonic array of values. Values are
 * assumed to increase or decrease monotonically, with no equal values.
 * This method is most efficient when called repeatedly for slightly 
 * changing search values; in such cases, the index returned from one 
 * call should be passed in the next.
 * <p>
 * Warning: this method does not ensure that the specified array is
 * monotonic; that check would be more costly than this search.
 * @param a the array of values, assumed to be monotonic.
 * @param x the value for which to search.
 * @param i the index at which to begin the search. If negative, this 
 *  method interprets this index as if returned from a previous call.
 * @return the index at which the specified value is found, or, if not
 *  found, -(i+1), where i equals the index at which the specified value 
 *  would be located if it was inserted into the monotonic array.
 */
int util::binarySearch(int* a, int alength, int x, int i) {
  int n = alength;
  int nm1 = n-1;
  int low = 0;
  int high = nm1;
  bool increasing = n<2 || a[0]<a[1];
  if (i<n) {
    high = (0<=i)?i:-(i+1);
    low = high-1;
    int step = 1;
    if (increasing) {
      for (; 0<low && x<a[low]; low-=step,step+=step)
        high = low;
      for (; high<nm1 && a[high]<x; high+=step,step+=step)
        low = high;
    } else {
      for (; 0<low && x>a[low]; low-=step,step+=step)
        high = low;
      for (; high<nm1 && a[high]>x; high+=step,step+=step)
        low = high;
    }
    if (low<0) low = 0;
    if (high>nm1) high = nm1;
  }
  if (increasing) {
    while (low<=high) {
      int mid = (low+high)>>1;
      int amid = a[mid];
      if (amid<x)
        low = mid+1;
      else if (amid>x)
        high = mid-1;
      else
        return mid;
    }
  } else {
    while (low<=high) {
      int mid = (low+high)>>1;
      int amid = a[mid];
      if (amid>x)
        low = mid+1;
      else if (amid<x)
        high = mid-1;
      else
        return mid;
    }
  }
  return -(low+1);
}

/**
 * Performs a binary search in a monotonic array of values. Values are
 * assumed to increase or decrease monotonically, with no equal values.
 * <p>
 * Warning: this method does not ensure that the specified array is
 * monotonic; that check would be more costly than this search.
 * @param a the array of values, assumed to be monotonic.
 * @param x the value for which to search.
 * @return the index at which the specified value is found, or, if not
 *  found, -(i+1), where i equals the index at which the specified value 
 *  would be located if it was inserted into the monotonic array.
 */
int util::binarySearch(int* a, int alength, int x) {
  return binarySearch(a,alength,x,alength);
}
  
////////////////////////////////////////////////////////////////////////////
// div
void util::div(float ra, float *ry, float *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = ra/ry[i];
}
void util::div(float ra, float **ry, float **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) rz[i2][i1] = ra/ry[i2][i1];
}
void util::div(float ra, float ***ry, float ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) rz[i3][i2][i1] = ra/ry[i3][i2][i1];
}

void util::mul(float ra, float *ry, float *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = ra*ry[i];
}

float* util::mul(float ra, float *rx, int n) {
  float* ry = zerofloat(n);
  mul(ra,rx,ry,n);
  return ry;
}

void util::mul(float ra, float **ry, float **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) mul(ra,ry[i2],rz[i2],n1);
}

////////////////////////////////////////////////////////////////////////////
// add
void util::add(float ra, float *ry, float *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = ra+ry[i];
}
void util::add(float ra, float **ry, float **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) rz[i2][i1] = ra+ry[i2][i1];
}
void util::add(float ra, float ***ry, float ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) rz[i3][i2][i1] = ra+ry[i3][i2][i1];
}
void util::add(float* rx, float *ry, float *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = rx[i]+ry[i];
}
void util::add(float** rx, float **ry, float **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) add(rx[i2],ry[i2],rz[i2],n1); 
}
void util::add(float*** rx, float ***ry, float ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) add(rx[i3],ry[i3],rz[i3],n1,n2); 
}
float* util::add(float* x, float* y, int n) {
  float* z = zerofloat(n);
  add(x,y,z,n);
  return z;
}
float** util::add(float** x, float** y, int n1, int n2) {
  float** z = zerofloat(n1,n2);
  add(x,y,z,n1,n2);
  return z;
}
float*** util::add(float*** x, float*** y, int n1, int n2, int n3) {
  float*** z = zerofloat(n1,n2,n3);
  add(x,y,z,n1,n2,n3);
  return z;
}
 
////////////////////////////////////////////////////////////////////////////
// sub
void util::sub(float ra, float *ry, float *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = ra-ry[i];
}
void util::sub(float ra, float **ry, float **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) rz[i2][i1] = ra-ry[i2][i1];
}
void util::sub(float ra, float ***ry, float ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) rz[i3][i2][i1] = ra-ry[i3][i2][i1];
}
float* util::sub(float ra, float* ry, int n) {
  float* rz = zerofloat(n);
  sub(ra,ry,rz,n);
  return rz;
}
float** util::sub(float ra, float** ry, int n1, int n2) {
  float** rz = zerofloat(n1,n2);
  sub(ra,ry,rz,n1,n2);
  return rz;
}
float*** util::sub(float ra, float*** ry, int n1, int n2, int n3) {
  float*** rz = zerofloat(n1,n2,n3);
  sub(ra,ry,rz,n1,n2,n3);
  return rz;
}
void util::sub(float* rx, float *ry, float *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = rx[i]-ry[i];
}
void util::sub(float** rx, float **ry, float **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) sub(rx[i2],ry[i2],rz[i2],n1); 
}
void util::sub(float*** rx, float ***ry, float ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) sub(rx[i3],ry[i3],rz[i3],n1,n2); 
}
float* util::sub(float* x, float* y, int n) {
  float* z = zerofloat(n);
  sub(x,y,z,n);
  return z;
}
float** util::sub(float** x, float** y, int n1, int n2) {
  float** z = zerofloat(n1,n2);
  sub(x,y,z,n1,n2);
  return z;
}
float*** util::sub(float*** x, float*** y, int n1, int n2, int n3) {
  float*** z = zerofloat(n1,n2,n3);
  sub(x,y,z,n1,n2,n3);
  return z;
}

void util::sub(double ra, double *ry, double *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = ra-ry[i];
}
void util::sub(double ra, double **ry, double **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) rz[i2][i1] = ra-ry[i2][i1];
}
void util::sub(double ra, double ***ry, double ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) 
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) rz[i3][i2][i1] = ra-ry[i3][i2][i1];
}
double* util::sub(double ra, double* ry, int n) {
  double* rz = zerodouble(n);
  sub(ra,ry,rz,n);
  return rz;
}
double** util::sub(double ra, double** ry, int n1, int n2) {
  double** rz = zerodouble(n1,n2);
  sub(ra,ry,rz,n1,n2);
  return rz;
}
double*** util::sub(double ra, double*** ry, int n1, int n2, int n3) {
  double*** rz = zerodouble(n1,n2,n3);
  sub(ra,ry,rz,n1,n2,n3);
  return rz;
}
void util::sub(double* rx, double *ry, double *rz, int n) {
  for (int i=0; i<n; ++i) rz[i] = rx[i]-ry[i];
}
void util::sub(double** rx, double **ry, double **rz, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) sub(rx[i2],ry[i2],rz[i2],n1); 
}
void util::sub(double*** rx, double ***ry, double ***rz, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) sub(rx[i3],ry[i3],rz[i3],n1,n2); 
}
double* util::sub(double* x, double* y, int n) {
  double* z = zerodouble(n);
  sub(x,y,z,n);
  return z;
}
double** util::sub(double** x, double** y, int n1, int n2) {
  double** z = zerodouble(n1,n2);
  sub(x,y,z,n1,n2);
  return z;
}
double*** util::sub(double*** x, double*** y, int n1, int n2, int n3) {
  double*** z = zerodouble(n1,n2,n3);
  sub(x,y,z,n1,n2,n3);
  return z;
}

////////////////////////////////////////////////////////////////////////////
// floats to doubles
double* util::float2double(float* f, int n) {
  double* g = new double[n];
  for (int i=0; i<n; ++i) g[i] = f[i];
  return g;
}
double** util::float2double(float** f, int n1, int n2) {
  double** g = new double*[n2];
  for (int i2=0; i2<n2; ++i2) 
    g[i2] = float2double(f[i2],n1);
  return g;
}
double*** util::float2double(float*** f, int n1, int n2, int n3) {
  double*** g = new double**[n3];
  for (int i3=0; i3<n3; ++i3) 
    g[i3] = float2double(f[i3],n1,n2);
  return g;
}

////////////////////////////////////////////////////////////////////////////
// doubles to floats
float* util::double2float(double* f, int n) {
  float* g = new float[n];
  for (int i=0; i<n; ++i) g[i] = f[i];
  return g;
}
float** util::double2float(double** f, int n1, int n2) {
  float** g = new float*[n2];
  for (int i2=0; i2<n2; ++i2) 
    g[i2] = double2float(f[i2],n1);
  return g;
}
float*** util::double2float(double*** f, int n1, int n2, int n3) {
  float*** g = new float**[n3];
  for (int i3=0; i3<n3; ++i3) 
    g[i3] = double2float(f[i3],n1,n2);
  return g;
}
  
////////////////////////////////////////////////////////////////////////////
// mulmat
// c[n][m] = a[n][l]*b[l][m]
void util::mulmat(float** a, float** b, float** c, int n, int m, int l) {
  int i, j, k;
  double** cc = new double*[n];
  for (i=0; i<n; ++i) cc[i] = new double[m];
  for (i=0; i<n; ++i)
    for (j=0; j<m; ++j) {
      cc[i][j] = 0;
      for (k=0; k<l; ++k) 
        cc[i][j] += a[i][k]*b[k][j];
    }
  for (i=0; i<n; ++i)
    for (j=0; j<m; ++j) c[i][j] = cc[i][j];
}
void util::mulmat(float** a, float *v, float *u, int n, int m) {
  double *uu = new double[m];
  int i;
  for (i=0; i<m; ++i) {
    uu[i] = 0;
    for (int j=0; j<n; ++j)
      uu[i] += a[i][j]*v[j];
  }
  for (i=0; i<m; ++i) u[i] = uu[i];
}

void util::mulmat(double **a, double **b, double **c, int n, int m, int l) {
  int i, j, k;
  for (i=0; i<n; ++i)
    for (j=0; j<m; ++j) {
      c[i][j] = 0;
      for (k=0; k<l; ++k) 
        c[i][j] += a[i][k]*b[k][j];
    }
}
void util::mulmat(double **a, double *v, double *u, int n, int m) {
  for (int i=0; i<m; ++i) {
    u[i] = 0;
    for (int j=0; j<n; ++j)
      u[i] += a[i][j]*v[j];
  }
}

///////////////////////////////////////////////////////////////////////////
// sum
float util::sum(float* x, int n) {
  float s=0;
  for (int i=0; i<n; ++i) s+=x[i];
  return s;
}
float util::sum(float** x, int n1, int n2) {
  float s=0;
  for (int i2=0; i2<n2; ++i2) s+=sum(x[i2],n1);
  return s;
}
float util::sum(float*** x, int n1, int n2, int n3) {
  float s=0;
  for (int i3=0; i3<n3; ++i3) s+=sum(x[i3],n1,n2);
  return s;
}
double util::sum(double* x, int n) {
  double s=0;
  for (int i=0; i<n; ++i) s+=x[i];
  return s;
}
double util::sum(double** x, int n1, int n2) {
  double s=0;
  for (int i2=0; i2<n2; ++i2) s+=sum(x[i2],n1);
  return s;
}
double util::sum(double*** x, int n1, int n2, int n3) {
  double s=0;
  for (int i3=0; i3<n3; ++i3) s+=sum(x[i3],n1,n2);
  return s;
}


///////////////////////////////////////////////////////////////////////////
// float <---> double

void util::float2double(float* x, double* y, int n1) {
  for (int i=0; i<n1; ++i) y[i] = x[i];
}

void util::float2double(float** x, double** y, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) float2double(x[i2],y[i2],n1);
}

void util::float2double(float*** x, double*** y, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) float2double(x[i3],y[i3],n1,n2);
}

void util::double2float(double* x, float* y, int n1) {
  for (int i=0; i<n1; ++i) y[i] = x[i];
}

void util::double2float(double** x, float** y, int n1, int n2) {
  for (int i2=0; i2<n2; ++i2) double2float(x[i2],y[i2],n1);
}

void util::double2float(double*** x, float*** y, int n1, int n2, int n3) {
  for (int i3=0; i3<n3; ++i3) double2float(x[i3],y[i3],n1,n2);
}

////////////////////////////////////////////////////////////////////////////
// correlation
double util::correlation(double** f, double** g, int n1, int n2) {
  double **ff = zerodouble(n1,n2);
  double **fg = zerodouble(n1,n2);
  double **gg = zerodouble(n1,n2);
  #pragma omp parallel for
  for (int i2=0; i2<n2; ++i2) 
    for (int i1=0; i1<n1; ++i1) {
      ff[i2][i1] = f[i2][i1]*f[i2][i1];
      fg[i2][i1] = f[i2][i1]*g[i2][i1];
      gg[i2][i1] = g[i2][i1]*g[i2][i1];
    }
  double fsum = sum(f,n1,n2);
  double gsum = sum(g,n1,n2);
  double fgsum = sum(fg,n1,n2);
  double ffsum = sum(ff,n1,n2);
  double ggsum = sum(gg,n1,n2);
  double nominator = n1*n2*fgsum-fsum*gsum;
  double denominator = (n1*n2*ffsum-fsum*fsum)*(n1*n2*ggsum-gsum*gsum);
  for (int i2=0; i2<n2; ++i2) {
    delete[] ff[i2]; 
    delete[] fg[i2];
    delete[] gg[i2];
  }
  delete [] ff;
  delete [] fg;
  delete [] gg;
  if (denominator>0) denominator = sqrt(denominator);
  else return 0;
  if (denominator!=0) return nominator/denominator;
  else return 0;
}

