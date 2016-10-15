#ifndef _UNITSPHERESAMPLING_H
#define _UNITSPHERESAMPLING_H

#include <iostream>

#include "ArrayMath.h"

using namespace std;

namespace util {
  /**
   * A roughly uniform sampling of the unit-sphere.
   * Maps integer sample indices i to points (x,y,z) on the unit sphere. The
   * sampling is only approximately uniform because sampled points on the
   * sphere are projections of points that uniformly sample an octahedron.
   * The corners of that octahedron are the points (1,0,0), (-1,0,0), (0,1,0),
   * (0,-1,0), (0,0,1), and (0,0,-1). Sampling near these corner points is
   * finer than sampling elsewhere on the unit sphere.
   * <p>
   * Positive sample indices correspond to points in the upper hemisphere
   * for which z&gt;=0. Negative sample indices correspond to points in the
   * lower hemisphere for which z&lt;=0. The sample index zero corresponds
   * to no point. Points on the equator (x,y,z=0) correspond to both negative
   * and positive sample indices.
   * <p>
   * Points with positive and negative indices are sampled symmetrically
   * about the center of the unit sphere. Specifically, if sample index i
   * corresponds to a point (x,y,z), then sample index -i corresponds to
   * the point (-x,-y,-z).
   *
   * @author Luming Liang, translated from Mines Java Toolkit
   * @version 2012.06.26
  */
  class UnitSphereSampling {
   
    ///////////////////////////////////////////////////////////////////////////
    // private

    // The unit sphere is projected onto an octahedron with corners that lie
    // on the sphere; that is, on the x, y, and z axes. The upper hemisphere
    // (z>=0) is projected onto the upper half of the octahedron.
    //
    // Imagine that this upper half is flattened into the x and y plane, and
    // let r and s the projected coordinates in this flat plane. The upper
    // hemisphere corresponds to a diamond in the r-s plane. The number of
    // points in this diamond is the number of points sampled for the upper
    // hemisphere, including the points on the equator for which z=0.
    //
    // The lower hemisphere (z<=0) is sampled in the same way, and its
    // mapping also includes points on the equator. The number of samples
    // in these diamonds (the sampling resolution) is limited by the number
    // of bits (nbits) in the signed integer sample indices.
    //
    // The r-s plane is sampled on an n by n grid, where n = 2*m+1 and m is
    // the number of sampling intervals along each of the positive r and s
    // axes. But not all samples in this grid are used.
    //
    // In the example below, for nbits = 6, m = 3, and n = 7, only the 25
    // samples marked with X correspond to sampled points. Samples marked
    // with 0 are unused, but are included in the grid for simplicity. The
    // indices for the upper hemisphere are in the range [1,25]. For the
    // lower hemisphere, indices are in the range [-25,-1]. The index 0 is
    // unused. This range of indices [-25,25] fits in a signed 6-bit integer
    // with range[-32:31].
    //
    //  s
    //  ^
    //  |
    //  3 |0 0 0 X 0 0 0
    //  2 |0 0 X X X 0 0
    //  1 |0 X X X X X 0
    //  0 |X X X X X X X ---> r
    // -1 |0 X X X X X 0
    // -2 |0 0 X X X 0 0
    // -3 |0 0 0 X 0 0 0
    // -------------------
    // -3 -2 -1 0 1 2 3
    //
    // Points on the equator correspond to the outermost points in the
    // the diamonds for the upper and lower hemispheres. These 4*m points
    // (for which z=0) appear in both tables. The number of unique points
    // sampled equals the number of points in the tables minus the 4*m
    // duplicate points.
    //
    // In a more practical example, nbits = 16, m = 127, and n = 255, with
    // sample indices in [-32513,-1] and [1,32513]. In this example, the
    // number of unique points sampled is 64518. This number is close to
    // but less than the maximum of 65536 points that could possibly be
    // represented in 16 bits.
    private: 
    int _m; // number of samples for positive r and s, not including zero
    int _n; // number of samples of r and s
    int _mindex; // maximum positive index
    int _nindex; // = _mindex+1; for negative indices
    int _npoint; // number of unique points
    double _d; // sampling interval for r and s = 1/m
    double _od; // one over sampling interval = m
    float** _pu; // table of points in upper hemisphere (z>=0)
    float** _pl; // table of points in lower hemisphere (z<=0)
    int** _ip; // table[n][n] of point indices
    
    static UnitSphereSampling* _uss16;
    static UnitSphereSampling* getUnitSphereSampling16() {
      if (_uss16==NULL) _uss16 = new UnitSphereSampling(16);
      return _uss16;
    }

    void initialize(int nbits) {
      //Check.argument(nbits>=4,"nbits>=4");
      //Check.argument(nbits<=32,"nbits<=32");

      // Sampling of the r-s plane with an n by n grid. Compute the
      // largest m such that the number of sample indices fits in a
      // signed integer with the specified number of bits. Note that
      // nbits-1 is the number of bits not counting the sign bit. The
      // upper limit on the largest positive index is 2^(nbits-1)-1.
      // The largest positive index is 1+2*m*(1+m), which also equals
      // the number (nindex) of positive indices.
      int indexLimit = (1<<(nbits-1))-1;
      int m = 1;
      while (1+2*m*(1+m)<=indexLimit)
        ++m;
      _m = --m;
      _n = 2*_m+1;

      // Sampling interval and its inverse in the r-s plane.
      _d = 1.0/_m;
      _od = _m;

      // Maximum positive index.
      _mindex = 1+2*_m*(1+_m);

      // Constant used for negative indices. Points in upper/lower
      // hemispheres are related by pu[index] = -pl[_nindex-index].
      _nindex = _mindex+1;

      // Number of unique sampled points. The number of points on the equator
      // is 4*m; these points (for which z=0) appear in the tables for both
      // upper and lower hemispheres. They have both positive and negative
      // indices, and are not counted twice here.
      _npoint = 2*_mindex-4*_m; // = 2+4*_m*_m

      //trace("m="+_m+" n="+_n+" mindex="+_mindex+" npoint="+_npoint);

      // Tables for points in upper and lower hemispheres.
      _pu = new float*[_nindex];
      _pl = new float*[_nindex];

      // Table of point indices.
      _ip = new int*[_n];
      for (int i=0; i<_n; ++i) _ip[i] = new int[_n];

      // For all sampled s on flattened octahedron, ...
      for (int is=0,js=-_m,index=0; is<_n; ++is,++js) {

        // Planar coordinate s and |s|.
        double s = js*_d;
        double as = (s>=0.0)?s:-s;

        // For all sampled r on flattened octahedron, ...
        for (int ir=0,jr=-_m; ir<_n; ++ir,++jr) {

          // Process only samples in the octahedral diamond corresponding
          // to the upper and lower hemispheres. Other points in the table
          // will be null.
          int jrs = abs(jr)+abs(js);
          if (jrs<=_m) {

            // Increment and store index in table.
            _ip[is][ir] = ++index;
            //trace("ir="+ir+" is="+is+" index="+index);

            // Planar coordinate r and |r|.
            double r = jr*_d;
            double ar = (r>=0.0)?r:-r;

            // Third coordinate t (t>=0) on octahedron.
            double t = max(0.0,1.0-ar-as);
            if (jrs==_m) t = 0.0;

            // Coordinates of point in upper hemisphere (z>=0).
            double scale = 1.0/sqrt(s*s+r*r+t*t);
            float x = (float)(r*scale);
            float y = (float)(s*scale);
            float z = (float)(t*scale);

            // Store coordinates in tables.
            float* pu = _pu[ index] = new float[3];
            float* pl = _pl[_nindex-index] = new float[3];
            pu[0] = x; pu[1] = y; pu[2] = z;
            pl[0] = -x; pl[1] = -y; pl[2] = -z;
          }
        }  
      }
    }

    float distanceOnSphere(float* p, float* q) {
      double x = p[0]+q[0];
      double y = p[1]+q[1];
      double z = p[2]+q[2];
      double d = x*x+y*y+z*z;
      if (d==0.0) {
        d = PI;
      } else if (d==4.0) {
        d = 0.0;
      } else {
        d = 2.0*atan(sqrt((4.0-d)/d));
      }
      return (float)d;
      //return (float)acos(p[0]*q[0]+p[1]*q[1]+p[2]*q[2]);
    }

    /**
     * Constructs a default sampling using sixteen bits per sample.
     * Errors for 16-bit samples are less than one degree of arc.
     * Where such errors are negligible, this default is both efficient
     * and convenient because points can be encoded as 16-bit shorts.
     */
    public:
    UnitSphereSampling() {
      initialize(16);
    }

    /**
     * Constructs a sampling for the specified number of bits.
     * Sample indices are signed integers with no more than this
     * number of bits, which includes the sign bit.
     * @param nbits the number of bits; 4 &lt;= nbits &lt;= 32 required.
     */
    UnitSphereSampling(int nbits) {
      initialize(nbits);
    }

    /**
     * Gets the number of points sampled on this unit sphere.
     * @return the number of points sampled.
     */
    int countSamples() {
      return _npoint;
    }
    
    /**
     * Gets the maximum sample index, a positive integer. The smallest
     * positive index is one. The smallest index is the negative of the
     * maximum index, and the largest negative index is minus one.
     * <p>
     * This number equals the number of points sampled in one hemisphere,
     * including points on the equator.
     * @return the maximum index.
     */
    int getMaxIndex() {
      return _mindex;
    }

    /**
     * Gets the sampled point for the specified index, which must be non-zero.
     * For efficiency, returns the array {x,y,z} of point coordinates
     * by reference, not by copy. These coordinates must not be modified.
     * @param index the index of the sampled point; must be non-zero.
     * @return array {x,y,z} of point coordinates; by reference, not by copy.
     */
    float* getPoint(int index) {
      return (index>=0)?_pu[index]:_pl[index+_nindex];
    }
    
    /**
     * Gets the index of the sampled point nearest to the specified point.
     * Here, the nearest sampled point is that nearest on the octahedron.
     * Returns a positive index for points in the upper hemisphere (z&gt;=0),
     * including points on the equator (z=0). Returns a negative index for
     * points in the lower hemisphere not on the equator (z&lt;0).
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point.
     * @param z z-coordinate of the point.
     * @return the sample index.
     */
    int getIndex(float x, float y, float z) {
      double ax = (x>=0.0f)?x:-x;
      double ay = (y>=0.0f)?y:-y;
      double az = (z>=0.0f)?z:-z;
      double scale = 1.0/(ax+ay+az);
      double r = x*scale;
      double s = y*scale;
      int ir = (int)(0.5+(r+1.0)*_od);
      int is = (int)(0.5+(s+1.0)*_od);
      int jr = ir-_m;
      int js = is-_m;
      if (jr+js>_m) {
        --ir;
        --is;
      } else if (-jr+js>_m) {
        ++ir;
        --is;
      } else if (-jr-js>_m) {
        ++ir;
        ++is;
      } else if (jr-js>_m) {
        --ir;
        ++is;
      }
      int index = _ip[is][ir];
      //assert index>0:"index>0";
      return (z>=0.0f)?index:index-_nindex;
    }

    /**
     * Gets the index of the sampled point nearest to the specified point.
     * Here, the nearest sampled point is that nearest on the octahedron.
     * @param xyz array {x,y,z} of point coordinates.
     * @return the sample index.
     */
    int getIndex(float* xyz) {
      return getIndex(xyz[0],xyz[1],xyz[2]);
    }

    /**
     * Gets an array {ia,ib,ic} of three sample indices for the spherical
     * triangle that contains the specified point. As viewed from outside
     * the sphere, the sampled points corresponding to the returned indices
     * are ordered counter-clockwise.
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point.
     * @param z z-coordinate of the point.
     * @return array of sample indices.
     */
    int* getTriangle(float x, float y, float z) {
      // Coordinates in r-s plane in [-1,1].
      double ax = (x>=0.0f)?x:-x;
      double ay = (y>=0.0f)?y:-y;
      double az = (z>=0.0f)?z:-z;
      double scale = 1.0/(ax+ay+az);
      double r = x*scale;
      double s = y*scale;
      // Integer grid indices in [0,2m]. These are the indices of the lower
      // left corner of the square in the grid that contains the point.
      double rn = (r+1.0)*_od;
      double sn = (s+1.0)*_od;
      int ir = (int)rn;
      int is = (int)sn;
      // Centered integer grid indices in [-m,m]. Useful for determining
      // which quadrant the point lies in, and whether indices lie outside
      // the sampled diamond. Points not outside satisfy |jr|+|js|<=m.
      int jr = ir-_m;
      int js = is-_m;
      // Adjust for points exactly on the equator in quadrant 1. The square
      // that contains the point must contain at least one triangle.
      if (jr+js==_m) {
        if (jr>0) {
          --jr;
          --ir;
        } else {
          --js;
          --is;
        }
      }
      // Fractional parts in [0,1). The grid square that contains the point
      // is split into two triangles. Squares in quadrants 1 and 3 have
      // lower-left and upper-right triangles. Squares in quadrants 2 and
      // 4 have upper-left and lower-right triangles. These fractional parts
      // are used to determine which triangle contains the point.
      double fr = rn-ir;
      double fs = sn-is;
      // Indices for sampled points of triangle.
      int ia,ib,ic;
      // If quadrant 1, ...
      if (jr>=0 && js>=0) {
        if (jr+js+2>_m || fr+fs<=1.0) {
          ia = _ip[is ][ir ]; // lower-left triangle
          ib = _ip[is ][ir+1];
          ic = _ip[is+1][ir ];
        } else {
          ia = _ip[is+1][ir+1]; // upper-right triangle
          ib = _ip[is+1][ir ];
          ic = _ip[is ][ir+1];
        }
      }
      // Else if quadrant 2, ...
      else if (jr<0 && js>=0) {
        if (-jr+js+1>_m || fr>=fs) {
          ia = _ip[is ][ir+1]; // lower-right triangle
          ib = _ip[is+1][ir+1];
          ic = _ip[is ][ir ];
        } else {
          ia = _ip[is+1][ir ]; // upper-left triangle
          ib = _ip[is ][ir ];
          ic = _ip[is+1][ir+1];
        }
      }
      // Else if quadrant 3, ...
      else if (jr<0 && js<0) {
        if (-jr-js>_m || fr+fs>=1.0) {
          ia = _ip[is+1][ir+1]; // upper-right triangle
          ib = _ip[is+1][ir ];
          ic = _ip[is ][ir+1];
        } else {
          ia = _ip[is ][ir ]; // lower-left triangle
          ib = _ip[is ][ir+1];
          ic = _ip[is+1][ir ];
        }
      }
      // Else if quadrant 4, ...
      else {
        if (jr+1-js>_m || fr<=fs) {
          ia = _ip[is+1][ir ]; // upper-left triangle
          ib = _ip[is ][ir ];
          ic = _ip[is+1][ir+1];
        } else {
          ia = _ip[is ][ir+1]; // lower-right triangle
          ib = _ip[is+1][ir+1];
          ic = _ip[is ][ir ];
        }
      }
      // All indices should be non-zero.
      if (ia==0 || ib==0 || ic==0) {
        cout<<"ia="<<ia<<" ib="<<ib<<" ic="<<ic<<endl;
        cout<<"x="<<x<<" y="<<y<<" z="<<z<<endl;
        cout<<"r="<<r<<" s="<<s<<endl;
        cout<<"ir="<<ir<<" is="<<is<<endl;
        cout<<"jr="<<jr<<" js="<<js<<endl;
        cout<<"fr="<<fr<<" fs="<<fs<<endl;
        cout<<"rn="<<rn<<" sn="<<sn<<endl;
       // assert false:"valid ia,ib,ic";
      }
      // Signs of indices depend on sign of z. Order the indices so that
      // points are in counter-clockwise order when viewed from outside
      // the sphere.
      int* ind = new int[3];
      if (z>=0.0f) {ind[0]=ia; ind[1]=ib; ind[2]=ic;}
      else {ind[0]=ia-_nindex; ind[1]=ic-_nindex; ind[2]=ib-_nindex;}
      return ind;
    }
    int* getTriangle(float* xyz) {
      return getTriangle(xyz[0],xyz[1],xyz[2]);
    }
    /**
     * Gets an array {wa,wb,wc} of three weights for a point in a spherical
     * triangle specified by sample indices of three points. The weights are
     * proportional to volumes of tetrahedra, and are used for interpolation.
     * Weights are non-negative and normalized so that their sum wa+wb+wc = 1.
     * <p>
     * For example, let p denote the specified point with coordinates {x,y,z},
     * and let o denote the center of the sphere with coordinates {0,0,0}.
     * Then the weight wa is proportional to the volume of the tetrahedron
     * formed by points p, b, c, and o.
     * @param x x-coordinate of the point.
     * @param y y-coordinate of the point.
     * @param z z-coordinate of the point.
     * @param iabc array {ia,ib,ic} of sample indices.
     * @return array {wa,wb,wc} of weights.
     */
    float* getWeights(float x, float y, float z, int* iabc) {
      float* pa = getPoint(iabc[0]);
      float* pb = getPoint(iabc[1]);
      float* pc = getPoint(iabc[2]);
      double xa = pa[0], ya = pa[1], za = pa[2];
      double xb = pb[0], yb = pb[1], zb = pb[2];
      double xc = pc[0], yc = pc[1], zc = pc[2];
      double wa = x*(yb*zc-yc*zb)+y*(zb*xc-zc*xb)+z*(xb*yc-xc*yb);
      double wb = x*(yc*za-ya*zc)+y*(zc*xa-za*xc)+z*(xc*ya-xa*yc);
      double wc = x*(ya*zb-yb*za)+y*(za*xb-zb*xa)+z*(xa*yb-xb*ya);
      if (wa<0.0) wa = 0.0;
      if (wb<0.0) wb = 0.0;
      if (wc<0.0) wc = 0.0;
      double ws = 1.0/(wa+wb+wc);
      float fa = (float)(wa*ws);
      float fb = (float)(wb*ws);
      float fc = (float)(wc*ws);
      float* w = new float[3];
      w[0] = fa; w[1] = fb; w[2] = fc;
      delete [] pa; delete [] pb; delete [] pc;
      return w;
    }
    float* getWeights(float* xyz, int* iabc) {
      return getWeights(xyz[0],xyz[1],xyz[2],iabc);
    }
    /**
     * Encodes specified points as 16-bit (short) indices.
     * @param x array of x-coordinates of points.
     * @param y array of y-coordinates of points.
     * @param z array of z-coordinates of points.
     * @return array of 16-bit (short) indices.
     */
    static short* encode16(float* x, float* y, float* z, int n1) {
      UnitSphereSampling* uss = getUnitSphereSampling16();
      short* s = new short[n1];
      for (int j=0; j<n1; ++j)
        s[j] = (short)uss->getIndex(x[j],y[j],z[j]);
      return s;
    }
    static short** encode16(float** x, float** y, float** z, int n1, int n2) {
      short** s = new short*[n2];
      for (int j=0; j<n2; ++j)
        s[j] = encode16(x[j],y[j],z[j],n1);
      return s;
    }
    static short*** encode16(float*** x, float*** y, float*** z, 
      int n1, int n2, int n3) {
      short*** s = new short**[n3];
      for (int j=0; j<n3; ++j)
        s[j] = encode16(x[j],y[j],z[j],n1,n2);
      return s;
    }
  };
}

#endif
