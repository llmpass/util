#ifndef VEC3D_H
#define VEC3D_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <string.h>
#include <stdarg.h>
#include <errno.h>

#include <iostream>
#include <fstream>
#include <vector>

#include "../lapack/tnt/tnt.h"

using namespace std;
using namespace TNT;

namespace util {
class Vec3d {
public:
  double x,y,z;

  Vec3d() {x=y=z=0;}
  Vec3d(double x,double y, double z) : x(x),y(y),z(z) {};
  Vec3d(double *f) : x(f[0]),y(f[1]),z(f[2]) {};
  
  void set(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  inline const double &operator[](const int i) const
  { return *(&x+i); }; 

  inline double &operator[](const int i)
  { return *(&x+i); }; 

  inline Vec3d &operator=(const Vec3d &b)
  { x = b.x; y = b.y; z = b.z; return *this;};
};

#define Epsilon 1E-5
#define Infinity HUGE_VAL

inline double Dot(const Vec3d &a, const Vec3d &b)
{ return a.x*b.x+a.y*b.y+a.z*b.z; };

inline Vec3d Product(const Vec3d &a, const Vec3d &b)
{ return Vec3d(a.x*b.x,a.y*b.y,a.z*b.z); };

inline Vec3d Cross(const Vec3d &a, const Vec3d &b){ 
  return Vec3d(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); 
};

inline Vec3d operator-(const Vec3d &v)
{ return Vec3d(-v.x,-v.y,-v.z); };

inline double Length(const Vec3d &v)
{ return sqrt(Dot(v,v)); };

inline Vec3d operator*(const double d, const Vec3d &v)
{ return Vec3d(d*v.x, d*v.y, d*v.z); };

inline Vec3d operator*(const Vec3d &v, const double d)
{ return Vec3d(d*v.x, d*v.y, d*v.z); };

inline void operator*=(Vec3d &v, const double d)
{ v.x*=d; v.y*=d; v.z*=d; };

inline void operator*=(Vec3d &v, const Vec3d &f)
{ v.x*=f.x; v.y*=f.y; v.z*=f.z; };

inline Vec3d operator/(const Vec3d &v, const double d)
{ return (1/d)*v; };

inline void operator/=(Vec3d &v, const double d)
{ v *= (1/d); };

inline bool operator==(Vec3d &v1, Vec3d &v2) 
{ return v1.x==v2.x && v1.y==v2.y && v1.z==v2.z;}

inline Vec3d operator+(const Vec3d &a, const Vec3d &b)
{ return Vec3d(a.x+b.x, a.y+b.y, a.z+b.z); };

inline Vec3d operator-(const Vec3d &a, const Vec3d &b)
{ return Vec3d(a.x-b.x, a.y-b.y, a.z-b.z); };

inline void Normalize(Vec3d &v)
{ if (Length(v)!=0) v *= (1.0f/Length(v)); };

inline Vector<double> ConjugateGradient(Sparse_Matrix<double>& A, Vector<double>
  b,Vector<double> x, int iMax=200, double eps=0.01) {
  int i = 0;
  Vector<double> r = b - A*x;
  Vector<double> d = r;
  double delta_new = dot_prod(r, r);
  double delta0 = delta_new, delta_old;
  Vector<double> q;
  double alpha, beta;
  Vector<double> xx(d.dim());
  xx = x;
  while (i < iMax && delta_new > eps*eps*delta0){
    q = A*d;
    alpha = delta_new / (dot_prod(d, q));
    xx = xx + alpha*d;
    if ((i%50)==0) r = b - A*xx;
    else r = r - alpha*q;
    delta_old = delta_new;
    delta_new = dot_prod(r, r);
    beta = delta_new / delta_old;
    d = r + beta*d;
    i++;
  }
  return xx;
}
}

#endif

