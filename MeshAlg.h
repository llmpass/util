#ifndef _MESHALG_H
#define _MESHALG_H

#include <queue>
#include <deque>
#include <time.h>
#include <iostream>
#include "Mesh.h"
#include "Andrzej/mesh.h"

using namespace std;

namespace util {
  class qitem {
    public:
    double d;
    int i, ancestor, gp;
    operator double() {return d;}
    qitem(double v, int ii, int a, int g) : d(v), i(ii), ancestor(a), gp(g) {}
  };

  void geoDis(Mesh& m, bool* mask, int* belongs, double*& d); 
  float linearMarching(Mesh& m, int i, double d, int i1, double d1, int i2);
  void mesh2Mesh(mesh& m, Mesh& M);
  void Mesh2mesh(Mesh& M, mesh& m);
  void disFace(Mesh&, int, int, double*&, int*);
  void disFace(Mesh&, int, int, double*&, vector<int>&);
  void segment(Mesh&, int, int, double*&);
  void avgGeoDis(Mesh& m);
  Mesh buildTorus(int, int);
  void gradient(Mesh& m, double** f, Vec3d*& g);
  void geoPath(Mesh& m, int s, int e, vector<int>& path);
  void findLastBif(vector<int>& path1, vector<int>& path2, int& b1, int& b2);
 
  /**
   * Compute the angle on the mesh from i0i1 to i0i2
   * the common vertex is i0
   */
  double accAng(Mesh& m, int i0, int i1, int i2) {
    int fInd, f1, f2;
    int i, j, k, ii, i3;
    double accumulateAng = 0, ang;
    while (i1!=i2) {
      m.findSharedFaces(i0,i1,f1,f2);
      // find the face that has orientation (i0-i1-i*)
      for (ii=0; ii<3; ++ii) {
        if (m.f(f1).v(ii)==i0 && m.f(f1).v((ii+1)%3)==i1) {fInd=f1; break;}
        if (m.f(f2).v(ii)==i0 && m.f(f2).v((ii+1)%3)==i1) {fInd=f2; break;}
      }
      // find the next vertex on this rotation
      for (ii=0; ii<3; ++ii) 
        if (m.f(fInd).v(ii)==i0 && m.f(fInd).v((ii+1)%3)==i1) 
          i3=m.f(fInd).v((ii+2)%3);
      // compute ang between i0i1 and i0i3
      double x01 = m.v(i1).x()-m.v(i0).x();
      double y01 = m.v(i1).y()-m.v(i0).y();
      double z01 = m.v(i1).z()-m.v(i0).z();
      Vec3d v01(x01,y01,z01);
      Normalize(v01);
      double x03 = m.v(i3).x()-m.v(i0).x();
      double y03 = m.v(i3).y()-m.v(i0).y();
      double z03 = m.v(i3).z()-m.v(i0).z();
      Vec3d v03(x03,y03,z03);
      Normalize(v03);
      accumulateAng += acos(Dot(v01,v03)); // [0,PI]
      //cout<<Dot(v01,v03)<<"  "<<acos(Dot(v01,v03))<<"  "<<accumulateAng<<endl;
      // move i1 to i3
      i1 = i3;
    }
    return accumulateAng;
  }

  double left(Mesh& m, int i0, int i1, int i2) {
    double accumulateAng = accAng(m,i0,i1,i2);
    if (accumulateAng<PI) return 1.0;
    else return -1.0;
  }
  
  /**
   * estimate gradient vector g of function f at every vertex of mesh m.
   * regression.
   */
  void gradient(Mesh& m, double* f, Vec3d*& g) {
    int ii, i, j, nv = m.NbVertex(), nf = m.NbFace();
    g = new Vec3d[nv];
    // gradient for each face
    Vec3d gf[nf];
    Vec3d vu, vv;
    double *u = zerodouble(3), *v = zerodouble(3);
    double *x = zerodouble(3), *y = zerodouble(3), *z = zerodouble(3);
    double *xx = zerodouble(3), *yy = zerodouble(3), *fv = zerodouble(3);
    double gx, gy;
    // for each face, build a 2d conformal triangle and local coordinates 
    for (i=0; i<nf; ++i) {
      for (j=0; j<3; ++j) {
        x[j] = m.v(m.f(i).v(j)).x();
        y[j] = m.v(m.f(i).v(j)).y();
        z[j] = m.v(m.f(i).v(j)).z();
        fv[j] = f[m.f(i).v(j)];
      }
      // para: 3 2d points xx[3], yy[3] <---> 3 3d points x[3], y[3], z[3]
      paraTri(x,y,z,xx,yy,vu,vv);
      // fit a linear expression in 2d domain to get the 2d gradient
      // solve gx*xx[j]+gy*yy[j]+c=fv[j] to obtain gx and gy. 
      if (xx[1]!=0) gx = (fv[1]-fv[0])/xx[1];
      else gx = 0;
      if (yy[2]!=0) gy = (fv[2]-fv[0]-gx*xx[2])/yy[2];
      else gy = 1;
      // convert (gx,gy) to 3d space
      gf[i] = vu*gx+vv*gy;
      Normalize(gf[i]);
    }
    delete [] u; delete [] v;
    delete [] x; delete [] y; delete [] z;
    delete [] xx; delete [] yy; 
    delete [] fv;
    for (ii=0; ii<nv; ++ii) {
      g[ii].set(0,0,0);
      for (i=0; i<m.v(ii).NbFaceNabors(); ++i) {
        // project the face gradient vector to the tangent plane of the vertex
        Point3 pi, pg;
        // naboring face index
        int fi = m.v(ii).fNabor(i);
        pi.set(m.v(ii).x(), m.v(ii).y(), m.v(ii).z());
        pg.set(m.v(ii).x()+gf[fi].x, m.v(ii).y()+gf[fi].y,
          m.v(ii).z()+gf[fi].z);
        Vec3d ni = m.v(ii).getNormal();
        Plane pl(pi,ni);
        Point3 pp = project(pg,pl);
        Vec3d vp = pp.toVec()-pi.toVec();
        Normalize(vp);
        g[ii] = g[ii]+vp;
      }
      Normalize(g[ii]);
    }
  }

  /**
   * Compute geodesic distance on a mesh from a set of vertices.
   * m is the mesh.
   * mask[i] = 1 indicates that ith vertex is a given vertex.
   */
  void geoDisFromVertexSet(Mesh& m, bool *mask, int* belongs, double*& d) {
    int n = m.NbVertex();
    for (int i=0; i<n; ++i) d[i] = 0;
    geoDis(m,mask,belongs,d);
  }

  /**
   * Compute geodesic distance on a mesh from a given vertex.
   * m is the mesh.
   * vertex v is a given vertex.
   */
  void geoDisFromVertex(Mesh& m, int v, double*& d) {
    int n = m.NbVertex();
    bool *mask = new bool[n];
    int* belongs = new int[n];
    for (int i=0; i<n; ++i) {
      d[i] = 0;
      mask[i] = false;
    }
    assert(v<n);
    mask[v] = true;
    belongs[v] = v;
    geoDisFromVertexSet(m,mask,belongs,d);
    delete [] mask;
    delete [] belongs;
  }
}

float util::linearMarching(Mesh& m, int i, double d, 
  int i1, double d1, int i2) {
  Vertex v  = m.v(i);
  Vertex v1 = m.v(i1);
  Vertex v2 = m.v(i2);
  // edge lengths of the triangle
  double a2 =  v.dis2(v1);
  double b2 =  v.dis2(v2);
  double c2 = v1.dis2(v2);
  double sina = (d1-d)/sqrt(a2); // -PI/2<alpha<PI/2
  if (sina>1) sina = 1;
  if (sina<-1) sina = -1;
  //float cosa = sqrt(1-a2); 
  double cosb = (a2+b2-c2)/(2*sqrt(a2)*sqrt(b2)); // 0<beta<PI
  double alpha = asin(sina);
  double beta  = acos(cosb);
  double d2 = d+sqrt(b2)*fabs(sin(alpha+beta));
  return d2;
}

void util::geoPath(Mesh& m, int s, int e, vector<int>& path) {
  double *d = zerodouble(m.NbVertex());
  geoDisFromVertex(m,s,d);
  vector<int> tempPath;
  // back tracking to find the path
  int i = e; 
  while (i!=s) {
    tempPath.push_back(i);
    int ii;
    // visit the naboring vertex of i, find the point with smallest distance
    double dis = 1000000; 
    for (int j=0; j<m.v(i).NbVertexNabors(); ++j) {
      int k = m.v(i).vNabor(j);
      if (d[k]<dis) {
        ii = k; 
        dis = d[k];
      }
    }
    i = ii;
  }
  tempPath.push_back(s);
  // reverse the path to ensure it begins from s, end at e
  path.resize(tempPath.size());
  for (i=0; i<tempPath.size(); ++i) path[i] = tempPath[tempPath.size()-1-i];
  if (s==2048 && (e==1603 || e==143)) 
    for (i=0; i<path.size(); ++i) cout<<path[i]<<"  ";
  cout<<endl;
}

/**
 * find the last bifurcation point of two pathes.
 * b1 stores the index of this point in path1
 * b2 stores the index of this point in path2
 **/
void util::findLastBif(vector<int>& path1, vector<int>& path2, 
  int& b1, int& b2) {
  int i1, i2, n1=path1.size(), n2=path2.size();
  i2 = n2-1;
  bool find = false;
  while (i2>=0) {
    for (i1=n1-1; i1>=0; --i1) { 
      if (path1[i1]==path2[i2]) {find = true; break;}
    }
    if (find) {b2 = i2; b1 = i1; break; }
    else i2--; 
  }  
  if (!find) {b1=-1; b2=-1;}
}

void util::geoDis(Mesh& m, bool* mask, int* belongs, double*& d) {
  priority_queue < qitem, deque<qitem>, greater<double> > pq;
  int i, j, k;
  for (i=0; i<m.NbVertex(); ++i)
    if (mask[i]) pq.push(qitem(0,i,0,belongs[i]));
  // start off by pushing the edge-propagated distance onto the priority queue
  // propagate the front to 1-ring nabors of marked vertices
  for (i=0; i<m.NbFace(); ++i)
    for (j=0; j<3; ++j) {
      int v1 = m.f(i).v((j+1)%3);
      int v2 = m.f(i).v((j+2)%3);
      if (mask[v1] && !mask[v2]) 
        pq.push(qitem(d[v1]+sqrt(m.v(v1).dis2(m.v(v2))),v2,0,belongs[v1]));
      if (!mask[v1] && mask[v2]) 
        pq.push(qitem(d[v1]+sqrt(m.v(v1).dis2(m.v(v2))),v1,0,belongs[v2]));
    }
  // main loop
  while (!pq.empty()) {
    // the current vertex is on the top of the priority queue
    qitem qi = pq.top();
    pq.pop();
    i = qi.i;
    if (mask[i]) continue;
    mask[i] = true;
    d[i] = qi.d;
    belongs[i] = qi.gp;
    Vertex& v = m.v(i);
    // update the queue by looping over the 1-ring nabors of the current
    // vertex
    for (j=0; j<v.NbVertexNabors(); ++j) {
      // nabor index
      int v1 = v.vNabor(j);
      if (!mask[v1]) {
        double dd = d[i]+sqrt(v.dis2(m.v(v1)));
        pq.push(qitem(dd,v1,0,belongs[i]));
      }
    }
    // update the queue by checking the 1-ring nabors of the current vertex,
    // if there is one nabor has been visited before, then we compute the
    // distance value of the remaining vertex in this triangle:
    // tri(v, v1, v2) 
    // v  = current vertex
    // v1 = marked nabor 
    // v2 = vertex with unkown distance value
    int nt = v.NbFaceNabors();
    // first find a nabor traingle of v that contains another marked nabor
    // vertex v1
    int vi1, vi2;
    for (j=0; j<nt; ++j) {
      int fi1 = v.fNabor(j);
      Face& f1 = m.f(fi1);
      // check if there is only one unknow vertex
      int numMarked = 0;
      for (k=0; k<3; ++k) {
        if (mask[f1.v(k)]) {
          // known vertex index
          if (f1.v(k)!=i) vi1 = f1.v(k);
          numMarked++;
        }
        // unknown vertex index
        else vi2 = f1.v(k);
      }
      if (numMarked!=2) continue;
      double dv2 = linearMarching(m,i,d[i],vi1,d[vi1],vi2);
      pq.push(qitem(dv2,vi2,0,belongs[vi1]));
    }
  }
}

void util::mesh2Mesh(mesh& m, Mesh& M) {
  int i, v1, v2, v3;
  int nf = m.faces(), nv = m.vertices();
  // dump faces
  for (i=0; i<nf; ++i) {
    const mesh_element *f = m.getface(i);
    // the vertices are indexed with odd numbers
    v1 = f->face[1]->ID;
    v2 = f->face[3]->ID;
    v3 = f->face[5]->ID;
    M.addFace(v1,v2,v3);
  }
  // dump vertices
  for (i=0; i<nv; ++i) {
    vec3dd v = m.vertex(i);
    M.addVertex(v[0],v[1],v[2]);
  }
  M.buildAdjacency();
  M.computeNormal();
}

void util::Mesh2mesh(Mesh& M, mesh& m) {
  int i, nf = M.NbFace(), nv = M.NbVertex();
  // dump faces
  for (i=0; i<nf; ++i) {
    Face& f = M.f(i);
    vector<int>* t = new vector<int>;
    t->push_back(f.v3());
    t->push_back(f.v1());
    t->push_back(f.v2());
    m.add_2Dmel(t);
  }
  // dump vertices
  delete [] m.v;
  m.v = new vec3dd[nv];
  for (i=0; i<nv; ++i) {
    m.v[i][0] = M.v(i).x();
    m.v[i][1] = M.v(i).y();
    m.v[i][2] = M.v(i).z();
  }
  m.finalize();
  m.compute_normals();
}

void util::disFace(Mesh& m, int seed, int patchNum, double*& geod, 
  int* maskInd) {
  int i, n = m.NbVertex(), ind, numMask = 0;
  geoDisFromVertex(m,seed,geod);
  maxv(geod,n,ind);
  //ind = seed;
  bool* mask = new bool[n];
  int* belongs = new int[n];
  maskInd[numMask] = ind;
  numMask++;
  for (i=0; i<n; ++i) mask[i] = false;
  mask[ind] = true;
  belongs[ind] = 0;
  geoDisFromVertexSet(m,mask,belongs,geod);
  while (numMask<patchNum) {
    maxv(geod,n,ind);
    maskInd[numMask] = ind;
    numMask++;
    for (i=0; i<n; ++i) mask[i] = false;
    for (i=0; i<numMask; ++i) {
      mask[maskInd[i]] = true;
      belongs[maskInd[i]] = i;
    }
    geoDisFromVertexSet(m,mask,belongs,geod);
  }
  delete [] mask;
  // set the group index for every vertex of m
  for (i=0; i<n; ++i) m.v(i).setGroup(belongs[i]);
  delete [] belongs;
}

/**
 * Segment mesh according to geodesic distances to seed vertices.
 **/
void util::disFace(Mesh& m, int seed, int patchNum, double*& geod, 
  vector<int>& maskInd) {
  int i, n = m.NbVertex(), ind;
  bool* mask = new bool[n];
  int* belongs = new int[n];
  maskInd.push_back(seed);
  for (i=0; i<n; ++i) mask[i] = false;
  mask[seed] = true;
  belongs[seed] = 0;
  geoDisFromVertexSet(m,mask,belongs,geod);
  while (maskInd.size()<patchNum) {
    maxv(geod,n,ind);
    maskInd.push_back(ind);
    for (i=0; i<n; ++i) mask[i] = false;
    for (i=0; i<maskInd.size(); ++i) {
      mask[maskInd[i]] = true;
      belongs[maskInd[i]] = i;
    }
    geoDisFromVertexSet(m,mask,belongs,geod);
  }
  delete [] mask;
  // set the group index for every vertex of m
  for (i=0; i<n; ++i) m.v(i).setGroup(belongs[i]);
  delete [] belongs;
}

// growing according to faces!!
void util::segment(Mesh& m, int seed, int patchNum, double*& geod) {
  // initialize complexs
  mesh me;
  Mesh2mesh(m,me);
  complexs* cp = new complexs(&me);
  int i, nf = m.NbFace(), v1, v2, v3, gn;
  vector<int> maskInd;
  disFace(m,seed,patchNum,geod,maskInd);
  // set the average geodesic distances of every cells
  // first 2d
  for (i=0; i<nf; ++i) cp->get_cell(i)->agd = m.f(i).getAgd();
  // then 1d
  for (i=0; i<me.edges(); ++i) {
    // 2 vertices
    v1 = me.getedge(i)->face[0]->ID;
    v2 = me.getedge(i)->face[1]->ID;
    cp->get_cell(i+me.faces())->agd = (m.v(v1).getAgd()+m.v(v2).getAgd())*0.5f;
  }
  // finally 0d
  for (i=0; i<m.NbVertex(); ++i) 
    cp->get_cell(i+me.edges()+me.faces())->agd = m.v(i).getAgd();
  // set the group index for every face of m by region growing
  bool* mask = new bool[nf];
  // compute currnt geodesic distances for every triangle faces
  double* fd = new double[nf];
  double onethree = 1.0/3;
  for (i=0; i<nf; ++i) {
    Face& f = m.f(i);
    fd[i] = (geod[f.v1()]+geod[f.v2()]+geod[f.v3()])*onethree;
    mask[i] = false;
  }
  priority_queue < qitem, deque<qitem>, greater<float> > pq;
  // seed vertices: point with geod = 0, stored in maskInd
  // starting triangles: tris contain seed vert
  for (i=0; i<maskInd.size(); ++i) {
    // find a tri that contains this vert
    int fi = m.v(maskInd[i]).fNabor(0);
    pq.push(qitem(fd[fi],fi,fi,i));
  }
  // main loop
  while (!pq.empty()) {
    qitem qi = pq.top();
    pq.pop();
    i = qi.i;
    if (mask[i]) continue;
    mask[i] = true;
    m.f(i).setGroup(qi.gp);
    // always merge the current cell to it's ancestors
    cell *ca = cp->get_cell(qi.ancestor), *ci = cp->get_cell(i);
    ca->merge(ci);
    // propogate the masked region from the current face to its unvisited
    // nabors. 
    for (int j=0; j<3; ++j) {
      int fk = m.f(i).getNabor(j);
      pq.push(qitem(fd[fk],fk,qi.ancestor,qi.gp));
    }
  }
  delete [] mask, fd;
  cp->remove_dead();
  // after merging, reset agd for 0d cells
  // To get a better correspondence between the 0d and 2d cells, for each 0d
  // cell, we set the average geodesic distance as the max agd of its
  // neighboring 2d cells
  int *m0, *i0;
  int n0 = cp->cells(0,m0,i0);
  for (i=0; i<n0; ++i) {
    cell* c0 = cp->get_cell(i0[i]);
    vector<int> c2Ar = c0->getNabor2dCellId();
    double sumAgd = 0; 
    for (int j=0; j<c2Ar.size(); ++j) 
      //if (cp->get_cell(c2Ar[j])->agd>maxAgd) 
        sumAgd += cp->get_cell(c2Ar[j])->agd;
    c0->agd = sumAgd/c2Ar.size();
  }
  m.setComplexs(cp);
  //cout<<*cp;
}

void util::avgGeoDis(Mesh& m) {
  int i, j, n = m.NbVertex(), sTime = time(NULL), tTime;
  int np = m.NbVertex()>200?200:m.NbVertex();
  // uniformly select np points to compute the geodesic distances
  Mesh mtemp = m;
  double* geod = new double[n];
  vector<int> vInd;
  disFace(mtemp,0,np,geod,vInd);
  double* agd = zerodouble(n);
  double* d = zerodouble(n);
  for (i=0; i<np; ++i) {
    geoDisFromVertex(m,vInd[i],d);
    for (j=0; j<n; ++j) agd[j] += d[j];
  }
  double np1 = 1.0f/np;
  for (i=0; i<n; ++i) m.v(i).setAgd(agd[i]*np1);
  for (i=0; i<m.NbFace(); ++i) {
    double s = 0;
    for (int j=0; j<3; ++j) 
      s += agd[m.f(i).v(j)];
    m.f(i).setAgd(s*0.3333f*np1);
  }
  delete [] d, agd, geod;  
  int eTime = time(NULL);
  cout<<"Runtime of computing agd is "<<eTime-sTime<<"s"<<endl;
}

using namespace util;
Mesh util::buildTorus(int n1, int n2) {
  Mesh mesh;
  int i1, i2, next1, next2, nV = n1*n2, NW, NE, SW, SE;
  // outer and inner radius
  float R=10.0, r=4.0;
  // evenly divide angles according to n1 and n2
  float v=PI*2.0/n1, u=PI*2.0/n2, rcosv, x, y, z;
  for (i1=0; i1<n1; ++i1) {
    z = r*sin(v*(i1+0.5)); rcosv = r*cos(v*(i1+0.5));
    next1 = (i1+1)%n1;
    for (i2=0; i2<n2; ++i2) {
      x = (R+rcosv)*cos(u*i2);
      y = (R+rcosv)*sin(u*i2);
      mesh.addVertex(x,y,z);
      next2 = (i2+1)%n2;
      // current quad
      // NW   NE
      // SW   SE
      NW = i1*n2+i2; NE = i1*n2+next2;
      SW = next1*n2+i2; SE = next1*n2+next2;
      // build 2 triangles for a quad (NW,SE,NE) and (NW,SW,SE) ccw
      mesh.addFace(NW,NE,SE);
      mesh.addFace(NW,SE,SW);
    }
  }
  Vec3d v3(0.9f,0.4f,0.2f);
  mesh.setColor(v3);
  mesh.setAlpha(0.7f);
  mesh.computeNormal();
  return mesh;
}

#endif
