#include "TriAlg.h"

/**
   * Parameterize a nondegenerate triangle.
   * Inputs: 3d coordiantes of triangle vertices: x, y, z arrays
   * Outputs: parameter u and v arrays; vu and vv 3d axis. 
   * Note: both u and v should be should be newed outside.
   */
void util::paraTri(double*& x, double*& y, double*& z, double*& u, double*& v,
  Vec3d& vu, Vec3d& vv) {
  Vec3d ve[3];
  for (int i=0; i<3; ++i) ve[i] = Vec3d(x[i],y[i],z[i]);
  Vec3d v01 = ve[1]-ve[0], v02 = ve[2]-ve[0];
  // normal
  Vec3d vn = Cross(v01,v02); Normalize(vn);
  // create u and v axis
  vu = v01; Normalize(vu);
  vv = Cross(vu,vn); Normalize(vv);
  // make sure vv adn v02 appear in the same half plane
  if (Dot(vv,v02)<0) vv = -vv;
  // get u and v arrays
  u[0] = 0; v[0] = 0;
  u[1] = Dot(v01,vu); v[1] = 0;
  u[2] = Dot(v02,vu); v[2] = Dot(v02,vv);
}

void util::paraTri(double*& x, double*& y, double*& z, double*& u, double*& v) {
  Vec3d vu, vv;
  paraTri(x,y,z,u,v,vu,vv);
}

/**
 * Area of a triangle.
 */
float util::area(float* x, float* y, float* z) {
  float l[3];
  for (int i=0; i<3; ++i)
    l[i] = sqrt((x[i]-x[(i+1)%3])*(x[i]-x[(i+1)%3])
               +(y[i]-y[(i+1)%3])*(y[i]-y[(i+1)%3])
               +(z[i]-z[(i+1)%3])*(z[i]-z[(i+1)%3]));
  float p = (l[0]+l[1]+l[2])*0.5f;
  return sqrt(p*(p-l[0])*(p-l[1])*(p-l[2]));
}

double util::area(double* x, double* y, double* z) {
  double l[3];
  for (int i=0; i<3; ++i)
    l[i] = sqrt((x[i]-x[(i+1)%3])*(x[i]-x[(i+1)%3])
               +(y[i]-y[(i+1)%3])*(y[i]-y[(i+1)%3])
               +(z[i]-z[(i+1)%3])*(z[i]-z[(i+1)%3]));
  double p = (l[0]+l[1]+l[2])*0.5f;
  return sqrt(p*(p-l[0])*(p-l[1])*(p-l[2]));
}

/**
 * Distortion metric between two triangles.
 * The distortion is defined as M = ||G-I||^2, where G = LL^T.
 *     |a b|   |ss.ss ss.st|
 * G = |b c| = |ss.st st.st|.
 * M = (a-1)^2+(c-1)^2+2b^2.
 * Note: the first triangle should NOT be degenerated.
 */
float util::distortion(double* x1, double* y1, double* z1, double* x2, 
  double* y2, double* z2) {
  double *u = new double[3], *v = new double[3];
  paraTri(x1,y1,z1,u,v);
  float area1 = 0.5f/area(x1,y1,z1);
  Vec3f q1(x2[0],y2[0],z2[0]);
  Vec3f q2(x2[1],y2[1],z2[1]);
  Vec3f q3(x2[2],y2[2],z2[2]);
  Vec3f ss = (q1*(v[1]-v[2])+q2*(v[2]-v[0])+q3*(v[0]-v[1]))*area1;
  Vec3f st = (q1*(u[1]-u[2])+q2*(u[2]-u[0])+q3*(u[0]-u[1]))*area1;
  delete [] u; delete [] v;
  float a = Dot(ss,ss), b = Dot(ss,st), c = Dot(st,st);
  float d;
  // d = a;
  d = sqrt(0.5f*(a+c)); d = fabs(d-1);
  //d = (a-1)*(a-1)+(c-1)*(c-1)+2*b*b;
  return d;
}

/**
 * Symmetric distortion metric.
 */
float util::distortSym(float* x1, float* y1, float* z1, float* x2, float* y2,
  float* z2) {
  // 1-->2: check if triangle 1 is degenerated or not
  int deg = degeneracy(x1,y1,z1);
  if (deg==2) {
    // perturb point
  }
  if (deg==1) {
    // perturb edge
  }
  // 2-->1
  return 0;
}

//float util::perturb
/**
 * Degeneracy of a triangle: 3: point, 2: edge, 1: triangle (nondegenerate).
 */
int util::degeneracy(float* x1, float* y1, float* z1) {
  int samePntNum = 1;
  for (int i=0; i<3; ++i) 
    if (x1[i]==x1[(i+1)%3] && y1[i]==y1[(i+1)%3] && z1[i]==z1[(i+1)%3]) 
      samePntNum++;
  return samePntNum;
}

