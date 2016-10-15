#ifndef _SPINIMAGE_H
#define _SPINIMAGE_H

#include "Mesh.h"
#include "ArrayMath.h"
#include "SegmentUnion.h"
#include "shot.h"
#include "../geometry/PointSet.h"
#include "../geometry/GeoAlg.h"
#include "../dsp/RecursiveGaussianFilter.h"

using namespace dsp;
using namespace geometry;

namespace util {
class SpinImage {
  public:
  int nImg, na, nb;
  double binSize;
  double ***im;

  /**
   * Constructor: we fix the image size, but leave the bin size free, so that
   * we could match spin images between objects with differet sizes.
   */
  SpinImage(int n1) {nImg = n1; na = n1; nb = 2*n1; im=NULL;}

  void compute(Mesh& m, int vi) {
    compute(m,vi,im[vi]);
  }

  void compute(Mesh& m, int vi, double**& img) {
    int nv = m.NbVertex();
    // the location and normal of the current vertex
    Vec3d p = m.v(vi).getV();
    Vec3d n = m.v(vi).getN(); 
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius()*0.8;
    // compute bin size: image width * binSize = support distance
    binSize = d/na*0.8;
    double binSize1 = 1.0f/binSize;
    // loop over all other vertices in mesh m to compute a and b
    for (int i=0; i<nv; ++i) {
      Vec3d x = m.v(i).getV();
      //if (Dot(x-p,x-p)>d*d) continue;
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double b = nxp;
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = b*binSize1+nImg;
      cout<<i<<"  "<<a<<"  "<<binSize<<endl;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      // bilinear interpolate to four nearest samples
      // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      cout<<ii<<"  "<<jj<<"  "<<na<<"  "<<nb<<endl;
      if (ii<na && jj<nb) img[ii  ][jj  ] += (1-s)*(1-t);
      if (ii+1<na && jj<nb) img[ii+1][jj  ] += s*(1-t);
      if (ii<na && jj+1<nb) img[ii  ][jj+1] += (1-s)*t;
      if (ii+1<na && jj+1<nb) img[ii+1][jj+1] += s*t; 
    }
  }

  void compute(Mesh& m, Vec3d& p, Vec3d& n, double**& img) {
    int nv = m.NbVertex();
    double d = m.bBox.radius();
    binSize = d/na;
    double binSize1 = 1.0f/binSize;
    // loop over all other vertices in mesh m to compute a and b
    for (int i=0; i<nv; ++i) {
      Vec3d x = m.v(i).getV();
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double b = nxp;
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = b*binSize1+nImg;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      // bilinear interpolate to four nearest samples
      // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      if (ii<na && jj<nb) img[ii  ][jj  ] += (1-s)*(1-t);
      if (ii+1<na && jj<nb) img[ii+1][jj  ] += s*(1-t);
      if (ii<na && jj+1<nb) img[ii  ][jj+1] += (1-s)*t;
      if (ii+1<na && jj+1<nb) img[ii+1][jj+1] += s*t; 
    }
    /*int n1 = nb, n2 = na;
    float** imgf = double2float(img,n1,n2);
    // difference of gaussian
    RecursiveGaussianFilter rgf1(1.0f);
    RecursiveGaussianFilter rgf2(3.0f);
    float **ss = zerofloat(n1,n2);
    float **ls = zerofloat(n1,n2);
    rgf1.apply00(imgf,ss,n1,n2);
    rgf2.apply00(imgf,ls,n1,n2);
    float **dif = zerofloat(n1,n2);
    sub(ss,ls,dif,n1,n2);
    float2double(dif,img,n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      delete [] imgf[i2]; 
      delete [] ss[i2]; delete [] ls[i2]; delete [] dif[i2]; 
    }
    delete [] imgf; 
    delete [] ss; delete [] ls; delete [] dif;*/ 
  }

  void compute(Mesh& m, Vec3d& p, Vec3d& n, double** img,
    double ratioSupport, double bs) {
    int nv = m.NbVertex();
    binSize = bs;
    double binSize1 = 1.0f/binSize;
    double d = m.bBox.radius();
    // loop over all other vertices in mesh m to compute a and b
    for (int i=0; i<nv; ++i) {
      Vec3d x = m.v(i).getV();
      double dis = Length(x-p);
      if (dis>d*ratioSupport) continue;
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double b = nxp;
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = b*binSize1+nImg;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      // bilinear interpolate to four nearest samples
      // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      if (ii<na && jj<nb) img[ii  ][jj  ] += (1-s)*(1-t);
      if (ii+1<na && jj<nb) img[ii+1][jj  ] += s*(1-t);
      if (ii<na && jj+1<nb) img[ii  ][jj+1] += (1-s)*t;
      if (ii+1<na && jj+1<nb) img[ii+1][jj+1] += s*t; 
    }
  }

  /**
   * Improved spin image: Zhang, Ong, and Foong, ICIP 2012
   * alpha and theta (angle)
   */
  void computeImprove(Mesh& m, Vec3d& p, Vec3d& n, double**& img, 
    double ratioSupport, double bs) {
    int nv = m.NbVertex();
    binSize = bs;
    double binSize1 = 1.0f/binSize;
    double d = m.bBox.radius();
    double angleSize = PI/na;
    double angleSize1 = 1.0/angleSize;
    double nx1, ny1, nz1;
    // loop over all other vertices in mesh m to compute a and b
    for (int i=0; i<nv; ++i) {
      Vec3d x = m.v(i).getV();
      Vec3d nx = m.v(i).getN();
      double dis = Length(x-p);
      if (dis>d*ratioSupport) continue;
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double D;
      if (nxp<Dot(nx,x-p)) D = 1;
      else D = -1;
      double dp = Dot(n,nx); 
      if (dp>1) dp = 1.0;
      if (dp<-1) dp = -1.0;
      double theta = D*acos(dp);
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = theta*angleSize1+nImg;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      // bilinear interpolate to four nearest samples
      // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      if (ii<na && jj<nb) img[ii  ][jj  ] += (1-s)*(1-t);
      if (ii+1<na && jj<nb) img[ii+1][jj  ] += s*(1-t);
      if (ii<na && jj+1<nb) img[ii  ][jj+1] += (1-s)*t;
      if (ii+1<na && jj+1<nb) img[ii+1][jj+1] += s*t; 
    }
  }

  void computeImprove(Mesh& m, Vec3d& p, Vec3d& n, double**& img) {
    int nv = m.NbVertex();
    double d = m.bBox.radius();
    binSize = d/na;
    double binSize1 = 1.0/binSize;
    double angleSize = PI/na;
    double angleSize1 = 1.0/angleSize;
    double nx1, ny1, nz1;
    // loop over all other vertices in mesh m to compute a and b
    for (int i=0; i<nv; ++i) {
      Vec3d x = m.v(i).getV();
      Vec3d nx = m.v(i).getN();
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double D;
      if (nxp<Dot(nx,x-p)) D = 1;
      else D = -1;
      double dp = Dot(n,nx); 
      if (dp>1) dp = 1.0;
      if (dp<-1) dp = -1.0;
      double theta = D*acos(dp);
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = theta*angleSize1+nImg;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      // bilinear interpolate to four nearest samples
      // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
      if (ii<na && jj<nb) img[ii  ][jj  ] += (1-s)*(1-t);
      if (ii+1<na && jj<nb) img[ii+1][jj  ] += s*(1-t);
      if (ii<na && jj+1<nb) img[ii  ][jj+1] += (1-s)*t;
      if (ii+1<na && jj+1<nb) img[ii+1][jj+1] += s*t; 
    }
  }

  void genStorage(Mesh& m) {
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    binSize = d/na;
    if (im==NULL) {
      im = new double**[m.NbVertex()];
      for (int i=0; i<m.NbVertex(); ++i) {
        im[i] = NULL;
        allocStorage(i);
      }
    }
  }

  void allocStorage(int vi) {
    if (im[vi]==NULL) {
      im[vi] = new double*[na];
      for (int i=0; i<na; ++i) {
        im[vi][i] = new double[nb];
        for (int j=0; j<nb; ++j) im[vi][i][j] = 0;
      }
    }
  }

  void genStorage(Mesh& m, int vi) {
    if (im==NULL) {
      double d = m.bBox.radius();
      binSize = d/na;
      im = new double**[m.NbVertex()];
      for (int i=0; i<m.NbVertex(); ++i) im[i] = NULL;
    }
    allocStorage(vi);
  }

  void releaseStorage(int vi) {
    if (im[vi]!=NULL) {
      for (int i=0; i<na; ++i) delete [] im[vi][i];
      delete [] im[vi];
      im[vi] = NULL;
    }
  }

  void releaseStorage(Mesh& m) {
    if (im!=NULL) {
      for (int i=0; i<m.NbVertex(); ++i) releaseStorage(i);
      delete [] im;
      im = NULL;
    }
  }

  /**
   * Generate the spin image of a traingular mesh.
   */
  void generate(Mesh& m) {
    //genStorage(m);
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    binSize = d/na;
    im = zerodouble(nb,na,m.NbVertex());
    // compute a and b for all vertices in the mesh
    for (int i=0; i<m.NbVertex(); ++i) compute(m,i);
  }

  /**
   * Compute the spin image of vertex vi.
   * Scan every triangle in the mesh, find out the range of a and b for each
   * triangle, then compute the overlapping area between the projected
   * triangle and the cylinder defined by a particular range of (a,b) pairs.
   */
  void computeScanTri(Mesh& m, int vi) {
    // the location and normal of the current vertex
    Vec3d p = m.v(vi).getV();
    Vec3d n = m.v(vi).getN();
    computeScanTri(m,p,n,im[vi]);
  }
  
  void computeScanTri(Mesh& m, int vi, double**& img) {
    // the location and normal of the current vertex
    Vec3d p = m.v(vi).getV();
    Vec3d n = m.v(vi).getN();
    computeScanTri(m,p,n,img);
  }

  void computeScanTri(Mesh& m, Vec3d& p, Vec3d& n, double**& img) {
    int i, nv = m.NbVertex(), nf = m.NbFace();
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    binSize = d/na;
    double binSize1 = 1.0f/binSize;
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = 0;
    Point3 tp[3]; // triangle points
    Vec3d x, vp1, vp2;
    // loop over triangles
    //#pragma omp parallel for private(tp,x,vp1,vp2)
    for (i=0; i<nf; ++i) {
      int iih = -1, jjh = -1; 
      int iil = na+1, jjl = nb+1;
      // compute the bin ranges associated to this triangle
      for (int j=0; j<3; ++j) {
        Vec3d x = m.v(m.f(i).v(j)).getV();
        tp[j] = Point3(x.x,x.y,x.z);
        double nxp = Dot(n,x-p);
        double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
        double b = nxp;
        double fi = a*binSize1;
        double fj = b*binSize1+nImg;
        int ii = (int)fi; 
        int jj = (int)fj;
        if (ii<iil) iil = ii;
        if (ii>iih) iih = ii;
        if (jj<jjl) jjl = jj;
        if (jj>jjh) jjh = jj;
      }
      // triangle polygon
      Polygon poly(tp,3);
      // loop over bins in [iil,iih]x[jjl,jjh]
      for (int ii=0; ii<=iih; ++ii)
        for (int jj=jjl; jj<=jjh; ++jj) {
          // for each bin, compute the prjoected contribution of the triangle.
          // first, create two planes parallel to the tangent plane of the
          // input vertex vi: n\dot(x-p) = jj-0.5 to jj+0.5
          vp1 = p+n*(jj-nImg)*binSize;
          vp2 = p+n*(jj-nImg+1)*binSize;
          Point3 po(p.x,p.y,p.z);
          Point3 p1(vp1.x,vp1.y,vp1.z);
          Point3 p2(vp2.x,vp2.y,vp2.z);
          Plane plt(po,n), pl1(p1,n), pl2(p2,n); // plt is the tangent plane
          // clip the polygon with these two planes
          //Polygon polyClip;
          //parallelClip(poly,pl1,pl2,polyClip);
          Polygon polyClip = parallelClip(poly,pl1,pl2);
          if (polyClip.isNull()) continue;
          // project the clipped polygon to the tangent plane
          Polygon polyProj = project(polyClip,plt);
          // second, generate two concentric circles on the tangent plane with
          // radii = ii-0.5f and ii+0.5f, the origin is at the input vertex
          double r1 = max(0.0,ii*1.0)*binSize;
          double r2 = (ii+1)*binSize;
          Circle c1(po,r1,plt), c2(po,r2,plt);
          // compute the overlapping area between the clipped and projected 
          // polygon and the ring bounded by two concentric circles.
          double area1 = areaOverCirclePolygon(c1,polyProj);
          double area2 = areaOverCirclePolygon(c2,polyProj);
          double overArea = area2-area1;
          // scale this area with the 1 over abs cosine of angle between 
          // the plane normal and the triangle normal
          Vec3d fn = m.f(i).getNormal();
          double factor = 1.0/fabs(Dot(n,fn));
          // finally, accumulate this area to the bin
          img[ii][jj] += overArea*factor;
        }
    }
    /*RecursiveGaussianFilter rgf(1.0f);
    float** imgf = double2float(img,nb,na);
    float **ss = zerofloat(nb,na);
    rgf.apply00(imgf,ss,nb,na);
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = ss[i][j];
    for (i=0; i<na; ++i) {
      delete [] ss[i]; delete [] imgf[i];
    }
    delete [] ss; delete [] imgf;*/
  }

  /////////////////////////////////////////////////////////////////////////
  // faster computation by fixing a better alpha range
  void computeScanTriFast(Mesh& m, int vi, double**& img) {
    // the location and normal of the current vertex
    Vec3d p = m.v(vi).getV();
    Vec3d n = m.v(vi).getN();
    computeScanTriFast(m,p,n,img);
  }

  void computeScanTriFast(Mesh& m, Vec3d& p, Vec3d& n, double**& img) {
    computeScanTriFast(m,p,n,img,1.0);
  }

  void computeScanTriFast(Mesh& m, Vec3d& p, Vec3d& n, double**& img, 
    double ratioSupport) {
    cout<<"in fast tri scanning"<<endl;
    int i, nv = m.NbVertex(), nf = m.NbFace();
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    binSize = d/na;
    double binSize1 = 1.0f/binSize;
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = 0;
    Point3 tp[3]; // triangle points
    Vec3d x, vp1, vp2;
    // segment union for computing contour in each horizontal line (beta)
    SegmentUnion su[nb];
    // loop over triangles
    //#pragma omp parallel for private(tp,x,vp1,vp2)
    for (i=0; i<nf; ++i) {
      Vec3d x0 = m.v(m.f(i).v(0)).getV();
      Vec3d x1 = m.v(m.f(i).v(1)).getV();
      Vec3d x2 = m.v(m.f(i).v(2)).getV();
      if (Length(x0-p)>d*ratioSupport && Length(x1-p)>d*ratioSupport
         && Length(x2-p)>d*ratioSupport) continue;
      int iih = -1, jjh = -1; 
      int iil = na+1, jjl = nb+1;
      // compute the beta ranges associated to this triangle
      for (int j=0; j<3; ++j) {
        Vec3d x = m.v(m.f(i).v(j)).getV();
        tp[j] = Point3(x.x,x.y,x.z);
        double nxp = Dot(n,x-p);
        double b = nxp;
        double fj = b*binSize1+nImg;
        int jj = (int)fj;
        if (jj<jjl) jjl = jj;
        if (jj>jjh) jjh = jj;
      }
      // triangle polygon
      Polygon poly(tp,3);
      // loop over beta \in [jjl, jjh] 
      for (int jj=jjl; jj<=jjh; ++jj) {
        // for each beta (jj), compute the range of alpha
        // first, create two planes parallel to the tangent plane of the
        // input vertex vi: n\dot(x-p) = jj to jj+1
        vp1 = p+n*(jj-nImg)*binSize;
        vp2 = p+n*(jj-nImg+1)*binSize;
        Point3 po(p.x,p.y,p.z);
        Point3 p1(vp1.x,vp1.y,vp1.z);
        Point3 p2(vp2.x,vp2.y,vp2.z);
        Plane plt(po,n), pl1(p1,n), pl2(p2,n); // plt is the tangent plane
        // clip the polygon with these two planes
        Polygon polyClip = parallelClip(poly,pl1,pl2);
        if (polyClip.isNull()) continue;
        // project the clipped polygon to the tangent plane
        Polygon polyProj = project(polyClip,plt);
        // compute the range of alpha (ii) 
        // first the lower bound: inside-->0, outside-->min dis to edges
        int iil, iih;
        int np = polyProj.n;
        double alpha=99999;
        bool isInside = inside(po,polyProj);
        if (isInside) iil=0;
        else {
          alpha=99999;
          // loop over segements
          for (int is=0; is<np; ++is) {
            double dist = dis(po,polyProj.s[is]);
            if (alpha>dist) alpha = dist;
          }
          iil = (int)(alpha*binSize1);
        }
        // then upper bound: max dis to vertices
        alpha=-1;
        for (int ip=0; ip<np; ++ip) {
          double dist = po.dis(polyProj.p[ip]);
          if (alpha<dist) alpha = dist;
        }
        iih = (int)(alpha*binSize1);
        for (int ii=iil; ii<=iih; ++ii) {
          // second, generate two concentric circles on the tangent plane with
          // radii = ii-0.5f and ii+0.5f, the origin is at the input vertex
          double r1 = max(0.0,ii*1.0)*binSize;
          double r2 = (ii+1)*binSize;
          Circle c1(po,r1,plt), c2(po,r2,plt);
          // compute the overlapping area between the clipped and projected 
          // polygon and the ring bounded by two concentric circles.
          double area1 = areaOverCirclePolygon(c1,polyProj);
          double area2 = areaOverCirclePolygon(c2,polyProj);
          double overArea = area2-area1;
          // scale this area with the 1 over abs cosine of angle between 
          // the plane normal and the triangle normal
          Vec3d fn = m.f(i).getNormal();
          double factor = 1.0/fabs(Dot(n,fn));
          // finally, accumulate this area to the bin
          img[ii][jj] = overArea*factor;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Compute the contour of the spin image by first computing the image
  void computeContour(Mesh& m, int vi, double**& img) {
    computeScanTriFast(m,vi,img);
    int n1 = nb, n2 = na;
    float** imgf = double2float(img,n1,n2); 
    float** silh = zerofloat(n1,n2);
    silhouette(imgf,silh,n1,n2);
    float2double(silh,img,n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      delete [] imgf[i2]; delete [] silh[i2];
    }
    delete [] imgf; delete [] silh;
  }

  void computeContour(Mesh& m, Vec3d& p, Vec3d& n, PointSet& ps) {
    int n1 = nb, n2 = na;
    double **img = zerodouble(n1,n2);
    computeScanTriFast(m,p,n,img);
    float** imgf = double2float(img,n1,n2); 
    float** silh = zerofloat(n1,n2);
    silhouette(imgf,silh,n1,n2);
    float2double(silh,img,n1,n2);
    vector<Point3> pset;
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1)
        if (silh[i2][i1]!=0) pset.push_back(Point3(i2,i1,0));
    for (int i2=0; i2<n2; ++i2) {
      delete [] imgf[i2]; delete [] silh[i2];
    }
    delete [] imgf; delete [] silh;
    Point3 pnts[pset.size()];
    for (int i=0; i<pset.size(); ++i) pnts[i] = pset[i];
    ps.set(pnts,pset.size());
  }
 
  void computeContourFast(Mesh& m, Vec3d& p, Vec3d& n, double**& img) {
    computeContourFast(m,p,n,img,1.0);
  }

  void computeContourFast(Mesh& m, Vec3d& p, Vec3d& n, double**& img, 
    double ratioSupport) {
    int i, nv = m.NbVertex(), nf = m.NbFace();
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    binSize = d/na;
    double binSize1 = 1.0f/binSize;
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = 0;
    Point3 tp[3]; // triangle points
    Vec3d x, vp1, vp2;
    // segment union for computing contour in each horizontal line (beta)
    SegmentUnion su[nb];
    // loop over triangles
    //#pragma omp parallel for private(tp,x,vp1,vp2)
    for (i=0; i<nf; ++i) {
      Vec3d x0 = m.v(m.f(i).v(0)).getV();
      Vec3d x1 = m.v(m.f(i).v(1)).getV();
      Vec3d x2 = m.v(m.f(i).v(2)).getV();
      if (Length(x0-p)>d*ratioSupport && Length(x1-p)>d*ratioSupport
         && Length(x2-p)>d*ratioSupport) continue;
      int iih = -1, jjh = -1; 
      int iil = na+1, jjl = nb+1;
      // compute the beta ranges associated to this triangle
      for (int j=0; j<3; ++j) {
        Vec3d x = m.v(m.f(i).v(j)).getV();
        tp[j] = Point3(x.x,x.y,x.z);
        double nxp = Dot(n,x-p);
        double b = nxp;
        double fj = b*binSize1+nImg;
        int jj = (int)fj;
        if (jj<jjl) jjl = jj;
        if (jj>jjh) jjh = jj;
      }
      // triangle polygon
      Polygon poly(tp,3);
      // loop over beta \in [jjl, jjh] 
      for (int jj=jjl; jj<=jjh; ++jj) {
        // for each beta (jj), compute the range of alpha
        // first, create two planes parallel to the tangent plane of the
        // input vertex vi: n\dot(x-p) = jj to jj+1
        vp1 = p+n*(jj-nImg)*binSize;
        vp2 = p+n*(jj-nImg+1)*binSize;
        Point3 po(p.x,p.y,p.z);
        Point3 p1(vp1.x,vp1.y,vp1.z);
        Point3 p2(vp2.x,vp2.y,vp2.z);
        Plane plt(po,n), pl1(p1,n), pl2(p2,n); // plt is the tangent plane
        // clip the polygon with these two planes
        Polygon polyClip = parallelClip(poly,pl1,pl2);
        if (polyClip.isNull()) continue;
        // project the clipped polygon to the tangent plane
        Polygon polyProj = project(polyClip,plt);
        // compute the range of alpha (ii) 
        // first the lower bound: inside-->0, outside-->min dis to edges
        int iil, iih;
        int np = polyProj.n;
        double alpha=99999;
        bool isInside = inside(po,polyProj);
        if (isInside) iil=0;
        else {
          alpha=99999;
          // loop over segements
          for (int is=0; is<np; ++is) {
            double dist = dis(po,polyProj.s[is]);
            if (alpha>dist) alpha = dist;
          }
          iil = (int)(alpha*binSize1);
        }
        // then upper bound: max dis to vertices
        alpha=-1;
        for (int ip=0; ip<np; ++ip) {
          double dist = po.dis(polyProj.p[ip]);
          if (alpha<dist) alpha = dist;
        }
        iih = (int)(alpha*binSize1);
        su[jj].addPair(iil,iih);
      }
    }
    // form contours by merging (union) segments
    int nseg, iil, iih;
    for (int jj=0; jj<nb; ++jj) {
      su[jj].unite();
      nseg = su[jj].segV.size();
      for (int iseg=0; iseg<nseg; ++iseg) {
        iil = (int)su[jj].segV[iseg].s;
        iih = (int)su[jj].segV[iseg].e;
        img[iil][jj] = 1;
        img[iih][jj] = 1;
      }
    }
  }
  
  void computeContourFast(Mesh& m, Vec3d& p, Vec3d& n, PointSet& ps) {
    computeContourFast(m,p,n,ps,1.0);
  }

  void computeContourFast(Mesh& m, Vec3d& p, Vec3d& n, PointSet& ps, 
    double ratioSupport) {
    int i, nv = m.NbVertex(), nf = m.NbFace();
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    binSize = d/na;
    double binSize1 = 1.0f/binSize;
    Point3 tp[3]; // triangle points
    Vec3d x, vp1, vp2;
    // segment union for computing contour in each horizontal line (beta)
    SegmentUnion su[nb];
    // loop over triangles
    //#pragma omp parallel for private(tp,x,vp1,vp2)
    for (i=0; i<nf; ++i) {
      Vec3d x0 = m.v(m.f(i).v(0)).getV();
      Vec3d x1 = m.v(m.f(i).v(1)).getV();
      Vec3d x2 = m.v(m.f(i).v(2)).getV();
      if (Length(x0-p)>d*ratioSupport && Length(x1-p)>d*ratioSupport
         && Length(x2-p)>d*ratioSupport) continue;
      int iih = -1, jjh = -1; 
      int iil = na+1, jjl = nb+1;
      // compute the beta ranges associated to this triangle
      for (int j=0; j<3; ++j) {
        Vec3d x = m.v(m.f(i).v(j)).getV();
        tp[j] = Point3(x.x,x.y,x.z);
        double nxp = Dot(n,x-p);
        double b = nxp;
        double fj = b*binSize1+nImg;
        int jj = (int)fj;
        if (jj<jjl) jjl = jj;
        if (jj>jjh) jjh = jj;
      }
      // triangle polygon
      Polygon poly(tp,3);
      // loop over beta \in [jjl, jjh] 
      for (int jj=jjl; jj<=jjh; ++jj) {
        // for each beta (jj), compute the range of alpha
        // first, create two planes parallel to the tangent plane of the
        // input vertex vi: n\dot(x-p) = jj to jj+1
        vp1 = p+n*(jj-nImg)*binSize;
        vp2 = p+n*(jj-nImg+1)*binSize;
        Point3 po(p.x,p.y,p.z);
        Point3 p1(vp1.x,vp1.y,vp1.z);
        Point3 p2(vp2.x,vp2.y,vp2.z);
        Plane plt(po,n), pl1(p1,n), pl2(p2,n); // plt is the tangent plane
        // clip the polygon with these two planes
        Polygon polyClip = parallelClip(poly,pl1,pl2);
        if (polyClip.isNull()) continue;
        // project the clipped polygon to the tangent plane
        Polygon polyProj = project(polyClip,plt);
        // compute the range of alpha (ii) 
        // first the lower bound: inside-->0, outside-->min dis to edges
        double iil, iih;
        int np = polyProj.n;
        double alpha=99999;
        bool isInside = inside(po,polyProj);
        if (isInside) iil=0;
        else {
          alpha=99999;
          // loop over segements
          for (int is=0; is<np; ++is) {
            double dist = dis(po,polyProj.s[is]);
            if (alpha>dist) alpha = dist;
          }
          iil = alpha;
        }
        // then upper bound: max dis to vertices
        alpha=-1;
        for (int ip=0; ip<np; ++ip) {
          double dist = po.dis(polyProj.p[ip]);
          if (alpha<dist) alpha = dist;
        }
        iih = alpha;
        su[jj].addPair(iil,iih);
      }
    }
    // form contours by merging (union) segments
    int nseg;
    double iil, iih, jjd;
    vector<Point3> pset;
    for (int jj=0; jj<nb; ++jj) {
      su[jj].unite();
      nseg = su[jj].segV.size();
      jjd = binSize*jj;
      for (int iseg=0; iseg<nseg; ++iseg) {
        iil = su[jj].segV[iseg].s;
        iih = su[jj].segV[iseg].e;
        pset.push_back(Point3(iil,jjd,0));
        pset.push_back(Point3(iih,jjd,0));
      }
    }
    // connect outer silheoutte
    double curS, curE, prevS, prevE, s, e; 
    bool isTop = true, isButtom = false;
    for (int jj=0; jj<nb; ++jj) {
      jjd = binSize*jj;
      int segSize = su[jj].segV.size();
      if (segSize!=0) {
        // need connect
        if (isTop) {
          s = su[jj].segV[0].s;
          e = su[jj].segV[0].e;
          rasterConnect(s,e,binSize,jjd,pset);
          prevS = su[jj].segV[0].s;
          prevE = su[jj].segV[segSize-1].e;
          isTop = false;
        }
        else {
          curS = su[jj].segV[0].s;
          curE = su[jj].segV[segSize-1].e;
          // connect curS-->prevS and curE-->prevE
          double sm = min(curS,prevS), la = max(curS,prevS);
          rasterConnect(sm,la,binSize,jjd,pset);
          sm = min(curE,prevE), la = max(curE,prevE);
          rasterConnect(sm,la,binSize,jjd,pset);
          prevS = su[jj].segV[0].s;
          prevE = su[jj].segV[segSize-1].e;
        }
      }
      else {
        if (jj>0 && su[jj-1].segV.size()>0) {
          // jj-1 is the buttom
          segSize = su[jj-1].segV.size();
          rasterConnect(su[jj-1].segV[0].s,su[jj-1].segV[segSize-1].e,
            binSize,(jj-1)*binSize,pset);
        }
        if (!isTop) break; 
      }
    }
    Point3 pnts[pset.size()];
    for (i=0; i<pset.size(); ++i) pnts[i] = pset[i];
    ps.set(pnts,pset.size());
  }

  // x1 must less or equal to x2
  void rasterConnect(double x1, double x2, double xInter, 
    double y, vector<Point3>& pset) {
    double x = x1+xInter;
    while (x<x2) {
      pset.push_back(Point3(x,y,0));
      x += xInter;
    }
  }

  void computeSegContour(Mesh& m, Vec3d& p, Vec3d& n, PointSet& ps) {
    double d = m.bBox.radius();
    computeSegContour(m,p,n,ps,1.0,d/na);
  }

  void computeSegContour(Mesh& m, Vec3d& p, Vec3d& n, PointSet& ps, 
    double ratioSupport, double bs) {
    int i, nv = m.NbVertex(), nf = m.NbFace();
    // compute support distance: the radius of the bounding box of an 
    double d = m.bBox.radius();
    // compute bin size: image width * binSize = support distance
    //binSize = d/na;
    binSize = bs;
    double binSize1 = 1.0f/binSize;
    Point3 tp[3]; // triangle points
    Vec3d x, vp1, vp2;
    // segment union for computing contour in each horizontal line (beta)
    SegmentUnion su[nb];
    // loop over triangles
    //#pragma omp parallel for private(tp,x,vp1,vp2)
    for (i=0; i<nf; ++i) {
      Vec3d x0 = m.v(m.f(i).v(0)).getV();
      Vec3d x1 = m.v(m.f(i).v(1)).getV();
      Vec3d x2 = m.v(m.f(i).v(2)).getV();
      if (Length(x0-p)>d*ratioSupport && Length(x1-p)>d*ratioSupport
         && Length(x2-p)>d*ratioSupport) continue;
      int iih = -1, jjh = -1; 
      int iil = na+1, jjl = nb+1;
      // compute the beta ranges associated to this triangle
      for (int j=0; j<3; ++j) {
        Vec3d x = m.v(m.f(i).v(j)).getV();
        tp[j] = Point3(x.x,x.y,x.z);
        double nxp = Dot(n,x-p);
        double b = nxp;
        double fj = b*binSize1+nImg;
        int jj = (int)fj;
        if (jj<jjl) jjl = jj;
        if (jj>jjh) jjh = jj;
      }
      // triangle polygon
      Polygon poly(tp,3);
      // loop over beta \in [jjl, jjh] 
      for (int jj=jjl; jj<=jjh; ++jj) {
        // for each beta (jj), compute the range of alpha
        // first, create a plane parallel to the tangent plane of the
        // input vertex vi: n\dot(x-p) = jj
        vp1 = p+n*(jj-nImg)*binSize;
        Point3 po(p.x,p.y,p.z);
        Point3 p1(vp1.x,vp1.y,vp1.z);
        Plane plt(po,n), pl1(p1,n); // plt is the tangent plane
        // intersect the polygon with this plane
        Segment s;
        bool isInter = intersect(pl1,poly,s);
        if (!isInter) continue;
        // project the segment to the tangent plane
        Segment sProj = project(s,plt);
        // compute the range of alpha (ii) 
        // first the lower bound: min dis to sProj
        double iil = dis(po,sProj);
        // then upper bound: max dis to vertices
        double dist1 = po.dis(sProj.s);
        double dist2 = po.dis(sProj.e);
        double iih = max(dist1,dist2);
        su[jj].addPair(iil,iih);
      }
    }
    // form contours by merging (union) segments
    int nseg;
    double iil, iih, jjd;
    vector<Point3> pset;
    for (int jj=0; jj<nb; ++jj) {
      su[jj].unite();
      nseg = su[jj].segV.size();
      jjd = binSize*jj;
      for (int iseg=0; iseg<nseg; ++iseg) {
        iil = su[jj].segV[iseg].s;
        iih = su[jj].segV[iseg].e;
        pset.push_back(Point3(iil,jjd,0));
        pset.push_back(Point3(iih,jjd,0));
      }
    }
    formContour(su,pset, 1);
    formContour(su,pset,-1);
    Point3 pnts[pset.size()];
    for (i=0; i<pset.size(); ++i) pnts[i] = pset[i];
    ps.set(pnts,pset.size());
  }

  // dir -- direction: 1--top to buttom -1--buttom to top
  void formContour(SegmentUnion* su, vector<Point3>& pset, int dir) {
    // loop over beta, connect intervals where it appears at current beta but
    // not appears at previous beta
    double curS, curE, prevS, prevE, s, e, jd; 
    bool isTop = true;
    vector<double> px;
    int start = dir>0 ? 0:nb-1;
    int end = dir>0 ? nb:-1;
    for (int j=start; j!=end; j+=dir) {
      jd = binSize*j;
      px.clear();
      int segSize = su[j].segV.size();
      if (segSize!=0) {
        // first raster connect every segments for this beta
        for (int i=0; i<segSize; ++i) {
          s = su[j].segV[i].s;
          e = su[j].segV[i].e;
          for (double ix=s; ix<e; ix+=binSize) px.push_back(ix);
        }
        // loop over every point in ps, decide which one should be included
        int np = px.size();
        for (int i=0; i<np; ++i) {
          if (isTop) {
            // every points are useful
            pset.push_back(Point3(px[i],jd,0));
            isTop = false;
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
            if (isNew) pset.push_back(Point3(px[i],jd,0));
          }
        }
      }
    }
  }

  void computeContour(Mesh& m, Vec3d& p, Vec3d& n, double**& img) {
    computeScanTri(m,p,n,img);
    int n1 = nb, n2 = na;
    float** imgf = double2float(img,n1,n2); 
    float** silh = zerofloat(n1,n2);
    silhouette(imgf,silh,n1,n2);
    float2double(silh,img,n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      delete [] imgf[i2]; delete [] silh[i2];
    }
    delete [] imgf; delete [] silh;
  }
  
  void silhouette(float** in, float** out, int n1, int n2) {
    zero(out,n1,n2);
    int i2l, i2h, i1l, i1h;
    for (int i2=0; i2<n2; ++i2) {
      i2l = max(0, i2-1); i2h = min(n2-1, i2+1);
      for (int i1=0; i1<n1; ++i1) {
        if (in[i2][i1]!=0) continue; 
        i1l = max(0, i1-1); i1h = min(n1-1, i1+1);
        // get pixels that have nonzero 8-neighbors 
        for (int ii2=i2l; ii2<=i2h; ++ii2)
          for (int ii1=i1l; ii1<=i1h; ++ii1) 
            if (in[ii2][ii1]!=0) out[ii2][ii1] = 1; //in[ii2][ii1];
      }
    }
  }

  /**
   * Generate the spin image of a triangular mesh by scan convert every
   * triangles in the mesh.
   */
  void generateRaster(Mesh& m) {
    genStorage(m);
    for (int i=0; i<m.NbVertex(); ++i) {
      computeScanTri(m,i);
      cout<<"generated "<<i<<"th spinimage"<<endl;
    }
  }

  /**
   * Approximate the continuous spin image using Monte Carlo sampling.
   * sv is the pre-generated samples.
   */
  void computeMonteCarlo(Mesh& m, Vec3d& p, Vec3d& n, 
    double**& img, double ratioSupport, double bs, vector<Vec3d>& sv) {
    int i, nv = m.NbVertex(), nf = m.NbFace();
    double d = m.bBox.radius();
    binSize = bs;
    double binSize1 = 1.0f/binSize;
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = 0;
     // loop over all samples to compute a and b
    for (int i=0; i<sv.size(); ++i) {
      Vec3d x = sv[i];
      double dis = Length(x-p);
      if (dis>d*ratioSupport) continue;
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double b = nxp;
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = b*binSize1+nImg;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      if (ii<na && jj<nb) img[ii  ][jj  ] += 1;//(1-s)*(1-t);
    }
  }

  void computeMonteCarlo(Mesh& m, int vi, double**& img, vector<Vec3d>& sv) {
    int i, nv = m.NbVertex(), nf = m.NbFace();
    // the location and normal of the current vertex
    Vec3d p = m.v(vi).getV();
    Vec3d n = m.v(vi).getN(); 
    double binSize1 = 1.0f/binSize;
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = 0;
     // loop over all samples to compute a and b
    for (int i=0; i<sv.size(); ++i) {
      Vec3d x = sv[i];
      double nxp = Dot(n,x-p);// n.(x-p)
      double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
      double b = nxp;
      // accumulate the contributions of this point to the image of vi: im[vi]
      double fi = a*binSize1, fj = b*binSize1+nImg;
      int ii = (int)fi, jj = (int)fj;
      double s = fi-ii, t = fj-jj;
      if (ii<na && jj<nb) im[vi][ii  ][jj  ] += 1;//(1-s)*(1-t);
    }
  }

  void computeMonteCarlo(Mesh& m, int vi, int no) {
    // generating samples
    vector<Vec3d> sv = m.genRandomSamples(no);
    computeMonteCarlo(m,vi,im[vi],sv);
  }

  void generateMonteCarlo(Mesh& m) {
    genStorage(m);
    int no = 10000;
    for (int i=0; i<m.NbVertex(); ++i) {
      cout<<"generating "<<i<<"th Monte Carlo spinimage"<<endl;
      computeMonteCarlo(m,i,no);
    }
  }

  void getSpinCoord(Vec3d n, Vec3d p, Vec3d x, double& fi, double& fj) {
    double binSize1 = 1.0f/binSize;
    double nxp = Dot(n,x-p);
    double a = sqrt(Dot(x-p,x-p)-nxp*nxp);
    double b = nxp;
    fi = a*binSize1;
    fj = b*binSize1+nImg;
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // Enhance the ridge of the spin image
  void computeRidge(Mesh& m, Vec3d& p, Vec3d& n, 
    double**& img, float inR, float outR) {
    computeScanTri(m,p,n,img);
    int n1 = nb, n2 = na;
    float** imgf = double2float(img,n1,n2);
    // difference of gaussian
    RecursiveGaussianFilter rgf1(inR);
    RecursiveGaussianFilter rgf2(outR);
    float **ss = zerofloat(n1,n2);
    float **ls = zerofloat(n1,n2);
    rgf1.apply00(imgf,ss,n1,n2);
    rgf2.apply00(imgf,ls,n1,n2);
    float **dif = zerofloat(n1,n2);
    sub(ss,ls,dif,n1,n2);
    // get ridge and valley
    /*int i2l, i2h, i1l, i1h;
    float **ridge = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      i2l = max(0, i2-1); i2h = min(n2-1, i2+1);
      for (int i1=0; i1<n1; ++i1) { 
        i1l = max(0, i1-1); i1h = min(n1-1, i1+1);
        // number of 8-neighbors that is smaller or larger than the pixel
        int sc=0, lc=0;
        for (int ii2=i2l; ii2<=i2h; ++ii2)
          for (int ii1=i1l; ii1<=i1h; ++ii1) {
            if (dif[ii2][ii1]<dif[i2][i1]) sc++;
            if (dif[ii2][ii1]>dif[i2][i1]) lc++;
          }
        if (sc>=6) ridge[i2][i1] = imgf[i2][i1]; 
      }
    }*/
    float2double(ss,img,n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      delete [] imgf[i2]; 
      delete [] ss[i2]; delete [] ls[i2]; delete [] dif[i2]; 
      //delete [] ridge[i2];
    }
    delete [] imgf; 
    delete [] ss; delete [] ls; delete [] dif; 
    //delete [] ridge;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Enhance the ridge of the spin image
  void computeRidge(Mesh& m, int vi, double**& img, float inR, float outR) {
    computeScanTri(m,vi,img);
    int n1 = nb, n2 = na;
    float** imgf = double2float(img,n1,n2);
    // difference of gaussian
    RecursiveGaussianFilter rgf1(inR);
    RecursiveGaussianFilter rgf2(outR);
    float **ss = zerofloat(n1,n2);
    float **ls = zerofloat(n1,n2);
    rgf1.apply00(imgf,ss,n1,n2);
    rgf2.apply00(imgf,ls,n1,n2);
    float **dif = zerofloat(n1,n2);
    sub(ss,ls,dif,n1,n2);
    // get ridge and valley
    /*int i2l, i2h, i1l, i1h;
    float **ridge = zerofloat(n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      i2l = max(0, i2-1); i2h = min(n2-1, i2+1);
      for (int i1=0; i1<n1; ++i1) { 
        i1l = max(0, i1-1); i1h = min(n1-1, i1+1);
        // number of 8-neighbors that is smaller or larger than the pixel
        int sc=0, lc=0;
        for (int ii2=i2l; ii2<=i2h; ++ii2)
          for (int ii1=i1l; ii1<=i1h; ++ii1) {
            if (dif[ii2][ii1]<dif[i2][i1]) sc++;
            if (dif[ii2][ii1]>dif[i2][i1]) lc++;
          }
        if (sc>=6) ridge[i2][i1] = imgf[i2][i1]; 
      }
    }*/
    float2double(dif,img,n1,n2);
    for (int i2=0; i2<n2; ++i2) {
      delete [] imgf[i2]; 
      delete [] ss[i2]; delete [] ls[i2]; delete [] dif[i2]; 
      //delete [] ridge[i2];
    }
    delete [] imgf; 
    delete [] ss; delete [] ls; delete [] dif; 
    //delete [] ridge;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Compute the partial ridges of a spin image by scanning a selected 
  // set of edges
  
  void generateScanEdge(Mesh& m) {
    genStorage(m);
    for (int i=0; i<m.NbVertex(); ++i) {
      cout<<"generating "<<i<<"th spinimage"<<endl;
      computeScanEdge(m,i,im[i]);
    }
  }
  
  void computeScanEdge(Mesh& m, int vi, double**& img) {
    int i, nv = m.NbVertex(), nf = m.NbFace(), ne = m.NbEdge();
    // the location and normal of the current vertex
    Vec3d p = m.v(vi).getV();
    Vec3d n = m.v(vi).getN(); 
    double binSize1 = 1.0f/binSize;
    for (i=0; i<na; ++i)
      for (int j=0; j<nb; ++j) img[i][j] = 0;
    int f1, f2;
    double ang1, ang2;
    // loop over edges
    for (i=0; i<ne; ++i) {
      // find two faces that share this edge
      m.findSharedFaces(m.e(i).getv1(),m.e(i).getv2(),f1,f2);
      // find two other (opposite) points in these two faces
      int v1 = m.findThirdPoint(m.e(i).getv1(),m.e(i).getv2(),f1);
      int v2 = m.findThirdPoint(m.e(i).getv1(),m.e(i).getv2(),f2);
      // compute the spin coordinates of these two points
      double iv1, jv1, iv2, jv2;
      Vec3d vx1 = m.v(v1).getV();
      Vec3d vx2 = m.v(v2).getV();
      getSpinCoord(n,p,vx1,iv1,jv1);
      getSpinCoord(n,p,vx2,iv2,jv2);
      // compute the bin ranges associated with this edge
      double ii1, jj1, ii2, jj2;
      Vec3d x1 = m.v(m.e(i).getv1()).getV();
      Vec3d x2 = m.v(m.e(i).getv2()).getV();
      getSpinCoord(n,p,x1,ii1,jj1);
      getSpinCoord(n,p,x2,ii2,jj2);
      // render the edge in the spin image only if two other points are located
      // at the same side of the edge (in the spin coordinates)
      Point3 e1(ii1,jj1,0);
      Point3 e2(ii2,jj2,0);
      Point3 e3(iv1,jv1,0);
      Point3 e4(iv2,jv2,0);
      // compute the angle between f1's normal and the current vertex normal
      double c1 = Dot(n,m.f(f1).getNormal());
      double c2 = Dot(n,m.f(f2).getNormal());
      ang1 = acos(Dot(n,m.f(f1).getNormal()))*180/PI;
      ang2 = acos(Dot(n,m.f(f2).getNormal()))*180/PI;
      // project the edge only if ang1<90 and ang2>90 or ang1>90 and ang1<90.
      //if (ang1<90 && ang2>90 || ang1>90 && ang2<90) {
      //if (fabs(ang1-ang2)>20) {
      if (sameSide(e1,e2,e3,e4)) {
        int iih = -1, jjh = -1; 
        int iil = na+1, jjl = nb+1;
        iil = (int)(ii1<ii2 ? ii1:ii2);
        iih = (int)(ii1<ii2 ? ii2:ii1);
        jjl = (int)(jj1<jj2 ? jj1:jj2);
        jjh = (int)(jj1<jj2 ? jj2:jj1);
        // edge segment
        Point3 ps(x1.x,x1.y,x1.z);
        Point3 pe(x2.x,x2.y,x2.z);
        Segment seg(ps,pe);
        //Polygon poly(tp,3);
        // loop over bins in [iil,iih]x[jjl,jjh]
        for (int ii=0; ii<=iih; ++ii)
          for (int jj=jjl; jj<=jjh; ++jj) {
            // for each bin, compute the length of this edge.
            // first, create two planes parallel to the tangent plane of the
            // input vertex vi: n\dot(x-p) = jj-0.5 to jj+0.5
            Vec3d vp1 = p+n*(jj-nImg)*binSize;
            Vec3d vp2 = p+n*(jj-nImg+1)*binSize;
            Point3 po(p.x,p.y,p.z);
            Point3 p1(vp1.x,vp1.y,vp1.z);
            Point3 p2(vp2.x,vp2.y,vp2.z);
            Plane plt(po,n), pl1(p1,n), pl2(p2,n); // plt is the tangent plane
            // clip the segment with these two planes
            Segment segClip = parallelClip(seg,pl1,pl2);
            // second, project the segment to the tangent plane
            Segment segProj = project(segClip,plt);
            // compute the partial length of the projected segment that lies 
            // in between the 2 circles
            // third, generate two concentric circles on the tangent plane with
            // radii = ii-0.5f and ii+0.5f, the origin is at the input vertex
            double r1 = max(0.0,ii*1.0)*binSize;
            double r2 = (ii+1)*binSize;
            Circle c1(po,r1,plt), c2(po,r2,plt);
            double l1 = lengthOverCircleSegment(c1,segProj); 
            double l2 = lengthOverCircleSegment(c2,segProj); 
            double lOver = l2-l1;
            // scale this length with the 1 over abs sine of angle between 
            // the plane normal and the segment vector
            Vec3d sv = seg.v;
            Normalize(sv);
            double cosa = Dot(n,sv);
            double sina = sqrt(1-cosa*cosa);
            double factor = 1.0/sina;
            // accumulate the scaled length to the bin
            img[ii][jj] += lOver*factor;
          }
      }
    } 
  }
};
}

#endif
