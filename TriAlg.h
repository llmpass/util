#ifndef _TRIALG_H
#define _TRIALG_H

#include "Vec3d.h"
#include "Vec3f.h"
using namespace std;

namespace util {
  /**
   * Parameterize a nondegenerate triangle.
   * Inputs: 3d coordiantes of triangle vertices: x, y, z arrays
   * Outputs: parameter u and v arrays; vu and vv 3d axis. 
   * Note: both u and v should be should be newed outside.
   */
  void paraTri(double*& x, double*& y, double*& z, double*& u, double*& v, 
    Vec3d& vu, Vec3d& vv);
  void paraTri(double*& x, double*& y, double*& z, double*& u, double*& v);
  /**
   * Area of a triangle.
   */
  float area(float* x, float* y, float* z);
  double area(double* x, double* y, double* z);
  /**
   * Distortion metric between two triangles.
   * The distortion is defined as M = ||G-I||^2, where G = LL^T.
   *     |a b|   |ss.ss ss.st|
   * G = |b c| = |ss.st st.st|.
   * M = (a-1)^2+(c-1)^2+2b^2.
   */
  float distortion(double* x1, double* y1, double* z1, double* x2, 
    double* y2, double* z2);
  float distortSym(float* x1, float* y1, float* z1, float* x2, float* y2,
    float* z2);
  int degeneracy(float* x1, float* y1, float* z1);
}

#endif
