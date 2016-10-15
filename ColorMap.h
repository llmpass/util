#ifndef _ColorMap_H
#define _ColorMap_H

#include "Vec3f.h"
#include "Vec3d.h"
#include "ArrayMath.h"

namespace util {
  Vec3f** jetColorMap(int n1, int n2, float** v, float min, float max) {
    Vec3f** cm = new Vec3f*[n2];
    float x, y, d = 1.0f/(max-min);
    for (int i2=0; i2<n2; ++i2) {
      cm[i2] = new Vec3f[n1];
      for (int i1=0; i1<n1; ++i1) {
        x = (v[i2][i1]-min)*d, y;
        if (x<0.125f) { // 0.0, 0.0, 0.5:1.0
          y = x/0.125f;
          cm[i2][i1] = Vec3f(0,0,0.5f+0.5f*y);
        } else if (x<0.375f) { // 0.0, 0.0:1.0, 1.0
          y = (x-0.125f)/0.25f;
          cm[i2][i1] = Vec3f(0,y,1.0f);
        } else if (x<0.625f) { // 0.0:1.0, 1.0, 1.0:0.0
          y = (x-0.375f)/0.25f;
          cm[i2][i1] = Vec3f(y,1.0f,1.0f-y);
        } else if (x<0.875f) { // 1.0, 1.0:0.0, 0.0
          y = (x-0.625f)/0.25f;
          cm[i2][i1] = Vec3f(1.0f,1.0f-y,0);
        } else { // 1.0:0.5, 0.0, 0.0
          y = (x-0.875f)/0.125f;
          cm[i2][i1] = Vec3f(1.0f-0.5f*y,0,0);
        }
      }
    }
    return cm;
  }
  
  Vec3f* jetColorMap(int n, float* v, float min, float max) {
    Vec3f* cm = new Vec3f[n];
    float x, y, d = 1.0f/(max-min);
    for (int i=0; i<n; ++i) {
      x = (v[i]-min)*d, y;
      if (x<0.125f) { // 0.0, 0.0, 0.5:1.0
        y = x/0.125f;
        cm[i] = Vec3f(0,0,0.5f+0.5f*y);
      } else if (x<0.375f) { // 0.0, 0.0:1.0, 1.0
        y = (x-0.125f)/0.25f;
        cm[i] = Vec3f(0,y,1.0f);
      } else if (x<0.625f) { // 0.0:1.0, 1.0, 1.0:0.0
        y = (x-0.375f)/0.25f;
        cm[i] = Vec3f(y,1.0f,1.0f-y);
      } else if (x<0.875f) { // 1.0, 1.0:0.0, 0.0
        y = (x-0.625f)/0.25f;
        cm[i] = Vec3f(1.0f,1.0f-y,0);
      } else { // 1.0:0.5, 0.0, 0.0
        y = (x-0.875f)/0.125f;
        cm[i] = Vec3f(1.0f-0.5f*y,0,0);
      }
    }
    return cm;
  }

  Vec3d** jetColorMap(int n1, int n2, double** v, double min, double max) {
    Vec3d** cm = new Vec3d*[n2];
    double x, y, d = 1.0f/(max-min);
    for (int i2=0; i2<n2; ++i2) {
      cm[i2] = new Vec3d[n1];
      for (int i1=0; i1<n1; ++i1) {
        x = (v[i2][i1]-min)*d, y;
        if (x<0.125f) { // 0.0, 0.0, 0.5:1.0
          y = x/0.125f;
          cm[i2][i1] = Vec3d(0,0,0.5f+0.5f*y);
        } else if (x<0.375f) { // 0.0, 0.0:1.0, 1.0
          y = (x-0.125f)/0.25f;
          cm[i2][i1] = Vec3d(0,y,1.0f);
        } else if (x<0.625f) { // 0.0:1.0, 1.0, 1.0:0.0
          y = (x-0.375f)/0.25f;
          cm[i2][i1] = Vec3d(y,1.0f,1.0f-y);
        } else if (x<0.875f) { // 1.0, 1.0:0.0, 0.0
          y = (x-0.625f)/0.25f;
          cm[i2][i1] = Vec3d(1.0f,1.0f-y,0);
        } else { // 1.0:0.5, 0.0, 0.0
          y = (x-0.875f)/0.125f;
          cm[i2][i1] = Vec3d(1.0f-0.5f*y,0,0);
        }
      }
    }
    return cm;
  }
  
  Vec3d* jetColorMap(int n, double* v, double min, double max) {
    Vec3d* cm = new Vec3d[n];
    double x, y, d = 1.0f/(max-min);
    for (int i=0; i<n; ++i) {
      x = (v[i]-min)*d, y;
      if (x<0.125f) { // 0.0, 0.0, 0.5:1.0
        y = x/0.125f;
        cm[i] = Vec3d(0,0,0.5f+0.5f*y);
      } else if (x<0.375f) { // 0.0, 0.0:1.0, 1.0
        y = (x-0.125f)/0.25f;
        cm[i] = Vec3d(0,y,1.0f);
      } else if (x<0.625f) { // 0.0:1.0, 1.0, 1.0:0.0
        y = (x-0.375f)/0.25f;
        cm[i] = Vec3d(y,1.0f,1.0f-y);
      } else if (x<0.875f) { // 1.0, 1.0:0.0, 0.0
        y = (x-0.625f)/0.25f;
        cm[i] = Vec3d(1.0f,1.0f-y,0);
      } else { // 1.0:0.5, 0.0, 0.0
        y = (x-0.875f)/0.125f;
        cm[i] = Vec3d(1.0f-0.5f*y,0,0);
      }
    }
    return cm;
  }

  Vec3d** grayColorMap(int n1, int n2, double** v, 
    double min, double max) {
    Vec3d** cm = new Vec3d*[n2];
    double x, d = 1.0f/(max-min);
    for (int i2=0; i2<n2; ++i2) {
      cm[i2] = new Vec3d[n1];
      for (int i1=0; i1<n1; ++i1) {
        x = (v[i2][i1]-min)*d;
        cm[i2][i1] = Vec3d(x,x,x);
      }
    }
    return cm;
  }
}

#endif
