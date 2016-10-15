#ifndef _GEOSPINIMAGE_H
#define _GEOSPINIMAGE_H

#include "Mesh.h"
#include "ArrayMath.h"
#include "SegmentUnion.h"
#include "../geometry/PointSet.h"
#include "../geometry/GeoAlg.h"
#include "../dsp/RecursiveGaussianFilter.h"

using namespace dsp;
using namespace geometry;

namespace util {

class GeoSpinImage {

public:
  int na, nb;
  double binSize1, binSize2;
  // array of geodesic distances 
  double** geoDis;
  double* avgGeoDis;
  double maxGeoDis;
  double maxAvgGeoDis, minAvgGeoDis;

  GeoSpinImage(int n1) {
    na = n1; nb = n1; 
    geoDis = NULL; avgGeoDis = NULL;
  }

  /**
   * Compute geodesic spin image using only vertex set.
   */
  void compute(Mesh& m, Vec3d& p, double**& img, double ratioSupport, 
    double binSize1, double support) {
    int i, j, nv = m.NbVertex(), nf = m.NbFace(), ii, jj;
    for (i=0; i<na; ++i) 
      for (j=0; j<nb; ++j) img[i][j] = 0;
    double binSize11 = 1.0/binSize1;
    // compute alpha bin size: image width * binSize = max avg geodesic 
    binSize2 = (maxAvgGeoDis-minAvgGeoDis)/na;
    double binSize21 = 1.0/binSize2;
    // find the nearest vertex of the given point
    int vp = m.nearestPnt(p.x, p.y, p.z);
    // compute the distance between the point and the nearest vertex
    double dp = (m.v(vp).x()-p.x)*(m.v(vp).x()-p.x);
    dp += (m.v(vp).y()-p.y)*(m.v(vp).y()-p.y);
    dp += (m.v(vp).z()-p.z)*(m.v(vp).z()-p.z);
    dp = sqrt(dp);
    // loop over all other vertices
    double a, b, fi, fj, s, t;
    for (i=0; i<nv; ++i) {
      a = avgGeoDis[i]-minAvgGeoDis;
      b = geoDis[vp][i]+dp;
      if (b>ratioSupport*support) continue;
      fi = a*binSize21; fj = b*binSize11;
      ii = (int)fi; jj = (int)fj;
      s = fi-ii; t = fj-jj;
      // bilinear interpolate to four nearest samples
      // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      if (ii<na && jj<nb) img[ii  ][jj  ] += (1-s)*(1-t);
      if (ii+1<na && jj<nb) img[ii+1][jj  ] += s*(1-t);
      if (ii<na && jj+1<nb) img[ii  ][jj+1] += (1-s)*t;
      if (ii+1<na && jj+1<nb) img[ii+1][jj+1] += s*t;
    }
  }
 
  /**
   * Approximate geodesic spin image using Monte Carlo sampling.
   */
  void computeMonteCarlo(Mesh& m, Vec3d& p, double**& img, double ratioSupport, 
    double binSize1, double support) {
  }
  
  void readGeoDis(double** gd, int nv) {
    int i, j;
    maxGeoDis = maxAvgGeoDis = -99999;
    minAvgGeoDis = 99999;
    geoDis = new double*[nv];
    avgGeoDis = zerodouble(nv);
    double nv1 = 1.0/nv;
    for (i=0; i<nv; ++i) {
      if (gd[i]!=NULL) {
        geoDis[i] = new double[nv];
        for (j=0; j<nv; ++j) {
          geoDis[i][j] = gd[i][j];
          avgGeoDis[i] += geoDis[i][j];
          if (maxGeoDis<geoDis[i][j]) maxGeoDis = geoDis[i][j];
        }
        avgGeoDis[i] *= nv1;
        if (avgGeoDis[i]>maxAvgGeoDis) maxAvgGeoDis = avgGeoDis[i];
        if (avgGeoDis[i]<minAvgGeoDis) minAvgGeoDis = avgGeoDis[i];
      }
    }
    cout<<"geodesic distances read"<<endl;
  }

  void readGeoDis(float** gd, int nv) {
    int i, j;
    maxGeoDis = maxAvgGeoDis = -99999;
    minAvgGeoDis = 99999;
    geoDis = new double*[nv];
    avgGeoDis = zerodouble(nv);
    float nv1 = 1.0/nv;
    for (i=0; i<nv; ++i) {
      if (gd[i]!=NULL) {
        geoDis[i] = new double[nv];
        for (j=0; j<nv; ++j) {
          geoDis[i][j] = gd[i][j];
          avgGeoDis[i] += geoDis[i][j];
          if (maxGeoDis<geoDis[i][j]) maxGeoDis = geoDis[i][j];
        }
        avgGeoDis[i] *= nv1;
        if (avgGeoDis[i]>maxAvgGeoDis) maxAvgGeoDis = avgGeoDis[i];
        if (avgGeoDis[i]<minAvgGeoDis) minAvgGeoDis = avgGeoDis[i];
      }
    }
    cout<<"geodesic distances read"<<endl;
  }

  /**
   * obtain avgDis only use the distance inside the ratio*support...
   **/
  void computeAvgGeoDis(Mesh& m, int nv, double ratio, double support) {
    int i, j, c;
    double d = support*ratio, dis, x, y, z;
    for (i=0; i<nv; ++i) {
      avgGeoDis[i] = 0;
      c = 0;
      for (j=0; j<nv; ++j) {
        x = m.v(i).x()-m.v(j).x();
        y = m.v(i).y()-m.v(j).y();
        z = m.v(i).z()-m.v(j).z();
        dis = x*x+y*y+z*z;
        if (geoDis[i][j]<d) {
          avgGeoDis[i] += geoDis[i][j];
          c++;
        }
      }
      avgGeoDis[i] /= c;
      if (avgGeoDis[i]>maxAvgGeoDis) maxAvgGeoDis = avgGeoDis[i];
      if (avgGeoDis[i]<minAvgGeoDis) minAvgGeoDis = avgGeoDis[i];
    }
    cout<<"average geodesic distances in a range computed"<<endl;
  }

  /**
   * geod[i][j] = geod[j][i];
   */
  void getGeoDis(Mesh& m) {
    int i, j, nv = m.NbVertex();
    maxGeoDis = maxAvgGeoDis = -999;
    minAvgGeoDis = 9999;
    geoDis = zerodouble(nv,nv);
    avgGeoDis = zerodouble(nv);
    double *temp = zerodouble(nv);
    double nv1 = 1.0/nv;
    for (i=0; i<nv; ++i) {
      geoDisFromVertex(m,i,temp);
      for (j=0; j<nv; ++j) {
        if (i<=j) geoDis[i][j] = temp[j];
        else geoDis[i][j] = (geoDis[j][i]+temp[j])*0.5;
        if (maxGeoDis<geoDis[i][j]) maxGeoDis = geoDis[i][j];
        avgGeoDis[i] += geoDis[i][j];
      }
      avgGeoDis[i] *= nv1;
      if (avgGeoDis[i]>maxAvgGeoDis) maxAvgGeoDis = avgGeoDis[i];
      if (avgGeoDis[i]<minAvgGeoDis) minAvgGeoDis = avgGeoDis[i];
    }
    cout<<"geodesic distances computed"<<endl;
    delete [] temp;
  }
};

}

#endif
