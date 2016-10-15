#ifndef _GEOSPINCONTOUR_H
#define _GEOSPINCONTOUR_H

#include "Mesh.h"
#include "ArrayMath.h"
#include "SegmentUnion.h"
#include "../geometry/PointSet.h"
#include "../geometry/PointSetf.h"
#include "../geometry/GeoAlg.h"
#include "../dsp/RecursiveGaussianFilter.h"

using namespace dsp;
using namespace geometry;

namespace util {

class GeoSpinContour {

public:
  int na, nb;
  float binSize1, binSize2;
  // array of geodesic distances 
  float** geoDis;
  float* avgGeoDis;
  float maxGeoDis;
  float maxAvgGeoDis, minAvgGeoDis;
  float maxA, minA;
  
  GeoSpinContour(int n1) {
    na = n1; nb = n1; 
    geoDis = NULL; avgGeoDis = NULL;
  }

  void compute(Mesh& m, Vec3d& p, PointSetf& ps) {
    compute(m, p, ps, 1);
  }

  void compute(Mesh& m, Vec3d& p, PointSetf& ps, double ratioSupport) {
    // compute support distance: greatest geodesic distance * 1.1
    double support = maxGeoDis*1.1;
    // compute beta bin size: image width * binSize = support distance
    binSize1 = support/na;
    compute(m,p,ps,ratioSupport,binSize1);
  }

  void compute(Mesh& m, Vec3d& p, PointSetf& ps, double ratioSupport, 
    double binSize1) {
    // compute support distance: greatest geodesic distance * 1.1
    double support = maxGeoDis*1.1;
    compute(m,p,ps,ratioSupport,binSize1,support);
  }

  void compute(Mesh& m, Vec3d& p, PointSetf& ps, double ratioSupport, 
    double binSize1, double support) {
    int i, j, nv = m.NbVertex(), nf = m.NbFace();
    double binSize11 = 1.0/binSize1;
    //getMinMax(m);
    // compute alpha bin size: image width * binSize = max avg geodesic 
    //binSize2 = (maxA-minA)/na;
    binSize2 = (maxAvgGeoDis-minAvgGeoDis)/na;
    double binSize21 = 1.0/binSize2;
    // find the nearest vertex of the given point
    int vp = m.nearestPnt(p.x, p.y, p.z);
    // compute the distance between the point and the nearest vertex
    double dp = (m.v(vp).x()-p.x)*(m.v(vp).x()-p.x);
    dp += (m.v(vp).y()-p.y)*(m.v(vp).y()-p.y);
    dp += (m.v(vp).z()-p.z)*(m.v(vp).z()-p.z);
    dp = sqrt(dp);
    // segment union for computing contour in each horizontal line (beta)
    SegmentUnion su[nb];
    // geodesic distances from the given point to triangle vertices
    double d[3];
    // average geodesic distances of the triangle vertics
    double ad[3];
    // three vertex indices
    int v[3];
    // status (inside or outside range) of thress vertices, true: in
    bool s[3];
    double lAgd1 = 1.0/(maxAvgGeoDis-minAvgGeoDis);
    double lGeoDis1 = 1.0/support;
    // loop over triangles
    //#pragma omp parallel for private(tp,x,vp1,vp2)
    for (i=0; i<nf; ++i) {
      // number of the vertex that are inside the range
      int numInside = 0; 
      for (j=0; j<3; ++j) {
        v[j] = m.f(i).v(j); 
        d[j] = geoDis[vp][v[j]]+dp;
        //ad[j] = avgGeoDis[v[j]]-avgGeoDis[vp]; // original: just use agd
        ad[j] = avgGeoDis[v[j]]; // original: just use agd
        s[j] = d[j]<support*ratioSupport;
        if (s[j]) numInside++; 
      }
      // if all outside, just discard
      if (numInside==0) continue;
      int iih = -1, jjh = -1; 
      int jjl = 9999, iil = 9999;
      // compute the beta ranges associated to this triangle
      // here, all beta are larger than 0!!!
      for (j=0; j<3; ++j) {
        int jj = d[j]*binSize11;
        if (jj<jjl) jjl = jj;
        if (jj>jjh) jjh = jj;
      }
      // loop over beta \in [jjl, jjh] 
      for (int jj=jjl; jj<=jjh; ++jj) {
        // for each beta (jj), compute the range of alpha
        // the alpha here is defined as the average geodesic distance
        // first, get the geodesic distance corresponding to jj
        double dj = jj*binSize1;
        if (dj>support) dj = support;
        // second, find which vertex of the give triangle has a larger
        // (smaller) geodesic distance than dj
        int smaller = 0, vd;
        for (j=0; j<3; ++j) {
          s[j] = d[j]<dj;
          if (s[j]) smaller++;
        }
        if (smaller==0 || smaller==3) continue; // should NOT happen...
        if (smaller==1) 
          for (j=0; j<3; ++j) if (s[j]) vd = j;
        if (smaller==2) 
          for (j=0; j<3; ++j) if (!s[j]) vd = j;
        // third, find two intersections (linear ratios) on edges
        // one on edge (vd, (vd+1)%3), another on edge (vd, (vd+2)%3)
        double r11 = (dj-d[(vd+1)%3])/(d[vd]-d[(vd+1)%3]);
        double r12 = 1-r11;
        double r21 = (dj-d[(vd+2)%3])/(d[vd]-d[(vd+2)%3]);
        double r22 = 1-r21;
        // linear interpolate average geodesic distances 
        double a1 = ad[vd]*r11+ad[(vd+1)%3]*r12;
        double a2 = ad[vd]*r21+ad[(vd+2)%3]*r22;
        // normalization
        a1 = (a1-minAvgGeoDis)*lAgd1;
        a2 = (a2-minAvgGeoDis)*lAgd1;
        if (fabs(p.x-0.062118)<1e-5 && fabs(p.y-0.042767)<1e-5) 
          cout<<a1<<"  "<<a2<<endl;
        su[jj].addPair(min(a1,a2),max(a1,a2));
        // avoid St9bad_alloc problem
        if (su[jj].nbSeg()>100) su[jj].unite();
      }
    }
    // form contours by merging (union) segments
    int nseg;
    double iil, iih, jjd;
    vector<Point3f> pset;
    for (int jj=0; jj<nb; ++jj) {
      su[jj].unite();
      nseg = su[jj].segV.size();
      jjd = binSize1*jj;
      // normalization
      jjd *= lGeoDis1;
      for (int iseg=0; iseg<nseg; ++iseg) {
        iil = su[jj].segV[iseg].s;
        iih = su[jj].segV[iseg].e;
        pset.push_back(Point3f(iil,jjd,0));
        pset.push_back(Point3f(iih,jjd,0));
      }
    }
    // scale bin size after normalization
    binSize1 = 1.0/nb; 
    binSize2 = 1.0/na;
    formContour(su,pset,-1,binSize1,binSize2);
    formContour(su,pset, 1,binSize1,binSize2);
    Point3f pnts[pset.size()];
    for (i=0; i<pset.size(); ++i) pnts[i].copy(pset[i]);
    ps.set(pnts,pset.size());
  }

  void releaseMem(Mesh& m) {
    int i, nv = m.NbVertex(); 
    /*if (geoDis!=NULL) {
      for (i=0; i<nv; ++i) 
        if (geoDis[i]!=NULL) delete [] geoDis[i];
      delete [] geoDis;
    }*/
    if (avgGeoDis!=NULL) delete [] avgGeoDis;
  }

  void readGeoDis(double** gd, int nv) {
    int i, j;
    maxGeoDis = maxAvgGeoDis = -99999;
    minAvgGeoDis = 99999;
    geoDis = new float*[nv];
    avgGeoDis = zerofloat(nv);
    float nv1 = 1.0/nv;
    for (i=0; i<nv; ++i) {
      if (gd[i]!=NULL) {
        geoDis[i] = new float[nv];
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
    geoDis = new float*[nv];
    avgGeoDis = zerofloat(nv);
    float nv1 = 1.0/nv;
    for (i=0; i<nv; ++i) {
      if (gd[i]!=NULL) {
        geoDis[i] = new float[nv];
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

  void useGeoDis(float** gd, int nv) {
    int i, j;
    maxGeoDis = maxAvgGeoDis = -99999;
    minAvgGeoDis = 99999;
    geoDis = gd;//new float*[nv];
    avgGeoDis = zerofloat(nv);
    float nv1 = 1.0/nv;
    for (i=0; i<nv; ++i) {
      for (j=0; j<nv; ++j) {
        avgGeoDis[i] += geoDis[i][j];
        if (maxGeoDis<geoDis[i][j]) maxGeoDis = geoDis[i][j];
      }
      avgGeoDis[i] *= nv1;
      if (avgGeoDis[i]>maxAvgGeoDis) maxAvgGeoDis = avgGeoDis[i];
      if (avgGeoDis[i]<minAvgGeoDis) minAvgGeoDis = avgGeoDis[i];
    }
    cout<<"geodesic distances read"<<endl;
  }

  /**
   * obtain avgDis only use the distance inside the ratio*support...
   **/
  void computeAvgGeoDis(Mesh& m, int nv, double ratio, double support) {
    int i, j, c;
    double d = support*ratio, dis, x, y, z;
    minAvgGeoDis = 99999; maxAvgGeoDis = -1;
    for (i=0; i<nv; ++i) {
      avgGeoDis[i] = 0;
      c = 0;
      for (j=0; j<nv; ++j) 
        if (geoDis[i][j]<d) {
          avgGeoDis[i] += geoDis[i][j];
          c++;
        }
      avgGeoDis[i] /= c;
      if (avgGeoDis[i]>maxAvgGeoDis) maxAvgGeoDis = avgGeoDis[i];
      if (avgGeoDis[i]<minAvgGeoDis) minAvgGeoDis = avgGeoDis[i];
    }
    cout<<"average geodesic distances in a range computed"<<endl;
  }

  void getMinMax(Mesh& m) {
    int i, j, nv = m.NbVertex();
    minA = maxf; maxA = -maxf;
    double a;
    for (i=0; i<nv; ++i) 
      for (j=0; j<nv; ++j) {
        a = avgGeoDis[i]/avgGeoDis[j];
        if (a<minA) minA = a;
        if (a>maxA) maxA = a;
      }
  }

  /**
   * geod[i][j] = geod[j][i];
   */
  void getGeoDis(Mesh& m) {
    int i, j, nv = m.NbVertex();
    maxGeoDis = maxAvgGeoDis = -999;
    minAvgGeoDis = 9999;
    geoDis = zerofloat(nv,nv);
    avgGeoDis = zerofloat(nv);
    double *temp = zerodouble(nv);
    float nv1 = 1.0/nv;
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

  // dir -- direction: 1--top to buttom -1--buttom to top
  void formContour(SegmentUnion* su, vector<Point3f>& pset, 
    int dir, double binSize1, double binSize2) {
    // loop over beta, connect intervals where it appears at current beta but
    // not appears at previous beta
    double curS, curE, prevS, prevE, s, e, jd; 
    bool isTop = true;
    vector<double> px; 
    int start = dir>0 ? 0:nb-1;
    int end = dir>0 ? nb:-1;
    for (int j=start; j!=end; j+=dir) {
      jd = binSize1*j;
      px.clear();
      int segSize = su[j].segV.size();
      if (segSize!=0) {
        // first raster connect every segments for this beta
        for (int i=0; i<segSize; ++i) {
          s = su[j].segV[i].s;
          e = su[j].segV[i].e;
          for (double ix=s; ix<e; ix+=binSize2) px.push_back(ix);
        }
        // loop over every point in ps, decide which one should be included
        int np = px.size();
        for (int i=0; i<np; ++i) {
          if (isTop) {
            // every points are useful
            pset.push_back(Point3f(px[i],jd,0));
            if (i==np-1) isTop = false;
          }
          else {
            // if this alpha does not lie in between any interval in the
            // previous beta, then we should include
            bool isNew = true;
            for (int i1=0; i1<su[j-dir].segV.size(); ++i1) {
              s = su[j-dir].segV[i1].s;
              e = su[j-dir].segV[i1].e;
              if (px[i]>s && px[i]<e) {
                isNew = false;
                break;
              }
            }
            if (isNew) pset.push_back(Point3f(px[i],jd,0));
          }
        }
      }
    }
  }
  
};

}

#endif
