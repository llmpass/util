#ifndef _MESH_H_
#define _MESH_H_

#include "Andrzej/complexs.h"
#include "../ann/include/ANN/ANN.h"

#include "Vec3d.h"
#include "Vec3f.h"
#include "vtkMath.h"
#include "../geometry/Point3.h"
#include "../geometry/GeoAlg.h"
#include "ArrayMath.h"
#include <vector>

using namespace geometry;
using namespace std;

namespace util {
class Edge {
  int v1, v2;
  public:
  Edge() {}
  Edge(int i1, int i2):v1(i1),v2(i2) {}
  void setv1(int i1) {v1=i1;}
  void setv2(int i2) {v2=i2;}
  int getv1() {return v1;}
  int getv2() {return v2;}
  bool isSame(Edge& e) {
    return (v1==e.v1 && v2==e.v2 || v1==e.v2 && v2==e.v1);
  }
};

class Face {
  int group;
  int _v1,_v2,_v3;
  int fNabors[3]; // fNabor[i] is opposite to v(i)
  float _nx,_ny,_nz;
  double _area;
  Vec3d _color;
  // average geodesic distance
  float agd;
  public: Face() {}
  Face(int v1, int v2,int v3) {
    _v1=v1; _v2=v2; _v3=v3;
  }
  int v1() {return _v1;}
  int v2() {return _v2;}
  int v3() {return _v3;}
  int v(int i) {
    if (i==0) return _v1;
    if (i==1) return _v2;
    if (i==2) return _v3;
  }
  float nx() {return _nx;}
  float ny() {return _ny;}
  float nz() {return _nz;}
  void setNormal(float& nx, float& ny, float& nz) {
    _nx=nx; _ny=ny; _nz=nz;
  }
  void setNormal(Vec3d& v) {
    _nx=v.x; _ny=v.y; _nz=v.z;
  }
  Vec3d getNormal() {
    return Vec3d(_nx,_ny,_nz);
  }
  void setColor(Vec3d& cv) {_color = cv;}
  void setColor(Vec3f& cv) {_color = Vec3d(cv.x,cv.y,cv.z);}
  Vec3d getColor() {return _color;}
  void setGroup(int g) {group = g;}
  int getGroup() {return group;}
  void setNabor(int i, int fi) {fNabors[i] = fi;}
  int getNabor(int i) {return fNabors[i];}
  bool contains(int v) {
    if (_v1==v || _v2==v || _v3==v) return true;
    return false;
  }
  void setAgd(float d) {agd = d;}
  float getAgd() {return agd;}
  void reorient() {
    int t = _v1; _v1 = _v2; _v2 = t;
    t = fNabors[0]; fNabors[0] = fNabors[1]; fNabors[1] = t;
  }
  void setArea(double area) {_area = area;}
  double getArea() {return _area;}
};

class Vertex {
  int group;
  float _x,_y,_z;
  float _nx,_ny,_nz;
  int _numFaceAround;
  Vec3d _color;
  float _alpha;
  // average geodesic distance
  float agd;
  // the edge nabors are edges starts at the current vertex
  vector<int> vNabors, fNabors, eNabors;
  public: 
  Vertex() {}  
  Vertex(float x, float y, float z) {
    _x=x; _y=y; _z=z; _numFaceAround=0;
    _nx=0; _ny=0; _nz=0;
    _color = Vec3d(1,0,0);
    _alpha = 0.7f;
  }
  Vec3d getV() {
    Vec3d v(_x,_y,_z);
    return v;
  }
  Vec3d getN() {
    Vec3d v(_nx,_ny,_nz);
    return v;
  }
  void getV(Vec3d& v) {
    v.set(_x,_y,_z);
  }
  void getN(Vec3d& n) {
    n.set(_nx,_ny,_nz);
  }
  void set(float x, float y, float z) {
    _x = x; _y = y; _z = z;
  }
  float x() {return _x;}
  float y() {return _y;}
  float z() {return _z;}
  float nx() {return _nx;}
  float ny() {return _ny;}
  float nz() {return _nz;}
  void addFaceNabor(int i) {fNabors.push_back(i);}
  void addVertexNabor(int i) {vNabors.push_back(i);}
  void addEdgeNabor(int i) {eNabors.push_back(i);}
  void clearFaceNabors() {fNabors.clear();}
  void clearVertexNabors() {vNabors.clear();}
  int NbFaceNabors() {return fNabors.size();}
  int NbVertexNabors() {return vNabors.size();}
  int fNabor(int i) {return fNabors[i];}
  int vNabor(int i) {return vNabors[i];}
  bool hasVertexNabor(int j) {
    for (int i=0; i<NbVertexNabors(); ++i) 
      if (vNabors[i]==j) return true;
    return false;
  }
  void setGroup(int g) {group = g;}
  int getGroup() {return group;}
  double dis2(Vertex& v) {
    double d2 = (_x-v.x())*(_x-v.x());
    d2 += (_y-v.y())*(_y-v.y());
    d2 += (_z-v.z())*(_z-v.z());
    return d2;
  }
  void setNormal(float& nx, float& ny, float& nz) {_nx=nx; _ny=ny; _nz=nz;}
  void setNormal(Vec3d& v) {_nx=v.x; _ny=v.y; _nz=v.z;}
  Vec3d getNormal() {return Vec3d(_nx,_ny,_nz);}
  void addNumFaceAroundByOne() {_numFaceAround++;}
  void addNormal(Vec3d& v) {_nx+=v.x; _ny+=v.y; _nz+=v.z;}
  void clearNormal() {_nx=0; _ny=0; _nz=0;}
  int getNumFaceAround() {return _numFaceAround;}
  void setColor(Vec3d& cv) {_color = cv;}
  void setColor(Vec3f& cv) {_color = Vec3d(cv.x,cv.y,cv.z);}
  Vec3d getColor() {return _color;}
  void setAlpha(float alpha) {_alpha = alpha;}
  float getAlpha() {return _alpha;}
  void setAgd(float d) {agd = d;}
  float getAgd() {return agd;}
  static Vertex& mid(Vertex v1, Vertex v2) {
    Vertex v((v1.x()+v2.x())*0.5f,(v1.y()+v2.y())*0.5f,(v1.z()+v2.z())*0.5f);
    return v;
  }
};

class BoundingBox {
  public:
  float xmin, xmax, ymin, ymax, zmin, zmax;
  float xsize, ysize, zsize, smax, smin;
  BoundingBox() {
    xmin = ymin = zmin = 99999;
    xmax = ymax = zmax = -99999;
  }
  void update(float x, float y, float z) {
    if (x<xmin) xmin = x; 
    if (x>xmax) xmax = x;
    if (y<ymin) ymin = y; 
    if (y>ymax) ymax = y; 
    if (z<zmin) zmin = z;
    if (z>zmax) zmax = z;
  }
  void computeSize() {
    xsize = xmax-xmin;
    ysize = ymax-ymin;
    zsize = zmax-zmin;
    smax = max(xsize,max(ysize,zsize));
    smin = min(xsize,min(ysize,zsize));
  }
  float radius() {
    float r = xsize*xsize+ysize*ysize+zsize*zsize;
    return sqrt(r);
  }
};

class Mesh {
  vector<Face> fTable;
  vector<Vertex> vTable;
  vector<Edge> eTable;
  // a pointer to a complex that records the segmentation
  complexs* cp;
  int *stack, sp, reorients;
  unsigned char* flags;
  // Ann kd tree 
  ANNkd_tree* pVertTree3;

  public: 
  BoundingBox bBox;
  Mesh() {cp = NULL; pVertTree3 = NULL;}
  ~Mesh() {
    if (cp) delete cp;
  }
  /**
   * Read off file.
   */
  void readOffFile(FILE *fp) {
    char off[3];
    fscanf(fp,"%s",off);
    int nv, nt, ng;
    fscanf(fp,"%d %d %d",&nv,&nt,&ng);
    fTable.resize(nt);
    vTable.resize(nv);
    int i, n, v1, v2, v3;
    float x, y, z;
    // read vertex table
    for (i=0; i<nv; ++i) {
      fscanf(fp,"%f %f %f",&x,&y,&z);
      vTable[i] = Vertex(x,y,z);
      bBox.update(x,y,z);
    }
    // read triangle table
    for (i=0; i<nt; ++i) {
      // shrec07, horse, mitAnim...
      fscanf(fp,"%d %d %d %d",&n,&v1,&v2,&v3);
      fTable[i] = Face(v1,v2,v3);
      // shrec2010
      //fscanf(fp,"%d %d %d",&v1,&v2,&v3);
      //fTable[i] = Face(v1-1,v2-1,v3-1);
    }
    buildAdjacency();
    computeNormal();
    bBox.computeSize();
  }
  void readOffFileWithColor(FILE *fp) {
    char coff[4];
    fscanf(fp,"%s",coff);
    int nv, nt, ng;
    fscanf(fp,"%d %d %d",&nv,&nt,&ng);
    fTable.resize(nt);
    vTable.resize(nv);
    int i, n, v1, v2, v3;
    float x, y, z;
    int cx, cy, cz, ca;
    // read vertex table
    for (i=0; i<nv; ++i) {
      fscanf(fp,"%f %f %f",&x,&y,&z);
      vTable[i] = Vertex(x,y,z);
      bBox.update(x,y,z);
      fscanf(fp,"%d %d %d %d",&cx,&cy,&cz,&ca);
      Vec3d v(cx*1.0/255.0,cy*1.0/255.0,cz*1.0/255.0);
      vTable[i].setColor(v);
    }
    // read triangle table
    for (i=0; i<nt; ++i) {
      // shrec07, horse, mitAnim...
      fscanf(fp,"%d %d %d %d",&n,&v1,&v2,&v3);
      fTable[i] = Face(v1,v2,v3);
      // shrec2010
      //fscanf(fp,"%d %d %d",&v1,&v2,&v3);
      //fTable[i] = Face(v1-1,v2-1,v3-1);
    }
    //buildAdjacency();
    //computeNormal();
    bBox.computeSize();
  }
  /**
    * Write off file.
    */
  void writeOffFile(string fn) {
    ofstream fs;
    fs.open(fn.c_str());
    fs<<"OFF"<<endl;
    fs<<vTable.size()<<" "<<fTable.size()<<" "<<0<<endl;
    // write vertices
    for (int i=0; i<NbVertex(); ++i) 
      fs<<v(i).x()<<" "<<v(i).y()<<" "<<v(i).z()<<endl;
    // write triangles
    for (int i=0; i<NbFace(); ++i) 
      fs<<"3  "<<f(i).v1()<<" "<<f(i).v2()<<" "<<f(i).v3()<<endl;
    fs.close();
  }
  /**
   * Read t file.
   */
  void readTFile(FILE *fp) {
    int nt, nv;
    fscanf(fp,"%d %d",&nt,&nv);
    fTable.resize(nt);
    vTable.resize(nv);
    // read triangle table
    int i, v1, v2, v3;
    for (i=0; i<nt; ++i) {
      fscanf(fp,"%d %d %d",&v1,&v2,&v3);
      fTable[i] = Face(v1,v2,v3);
    }
    // read vertex table
    float x, y, z;
    for (i=0; i<nv; ++i) {
      fscanf(fp,"%f %f %f",&x,&y,&z);
      vTable[i] = Vertex(x,y,z);
      bBox.update(x,y,z);
    }
    buildAdjacency();
    computeNormal();
    bBox.computeSize();
  }
  /**
    * Write t file.
    */
  void writeTFile(string fn) {
    ofstream fs;
    fs.open(fn.c_str());
    fs<<fTable.size()<<" "<<vTable.size()<<endl<<endl;
    // write triangles
    for (int i=0; i<NbFace(); ++i) 
      fs<<f(i).v1()<<" "<<f(i).v2()<<" "<<f(i).v3()<<endl;
    // write vertices
    for (int i=0; i<NbVertex(); ++i) 
      fs<<v(i).x()<<" "<<v(i).y()<<" "<<v(i).z()<<endl;
    fs.close();
  }
  complexs* getComplexs() {return cp;}
  void setComplexs(complexs* c) {cp = c;}

  /**
   * the center of mass of a mesh.
   */
  void center(Point3& p) {
    float x=0, y=0, z=0;
    int n = NbVertex();
    for (int i=0; i<n; ++i) {
      x += v(i).x();
      y += v(i).y();
      z += v(i).z();
    }
    float n1 = 1.0f/n;
    p.set(x*n1,y*n1,z*n1);
  }

  void addFace(int v1, int v2, int v3) {fTable.push_back(Face(v1,v2,v3));}
  void addVertex(float x, float y, float z) {vTable.push_back(Vertex(x,y,z));}
  int NbFace() {return fTable.size();}
  int NbEdge() {return eTable.size();}
  int NbVertex() {return vTable.size();}
  Face& f(int i) {return fTable[i];} 
  Vertex& v(int i) {return vTable[i];}
  Edge& e(int i) {return eTable[i];}
  void setColor(Vec3d& cv) {
    for (int i=0; i<NbVertex(); ++i)
      vTable[i].setColor(cv);
  }
  void setColor(Vec3f& cv) {
    for (int i=0; i<NbVertex(); ++i)
      vTable[i].setColor(cv);
  }
  Vec3d getColor() {return vTable[0].getColor();}
  void setAlpha(float alpha) {
    for (int i=0; i<NbVertex(); ++i) 
      vTable[i].setAlpha(alpha);
  }
  float getAlpha() {return vTable[0].getAlpha();}

  ANNkd_tree* buildAnnKdTree() {
    int n = vTable.size();
    ANNpointArray vArray = annAllocPts(n,3);
    for (int i=0; i<n; ++i) {
      vArray[i][0] = vTable[i].x();
      vArray[i][1] = vTable[i].y();
      vArray[i][2] = vTable[i].z();
    }
    //Build kd Tree for nearest distance searching
    return new ANNkd_tree(vArray,n,3);
  }

  /////////////////////////////////////////////////////////////////////////
  // nearest neighbor (vertices) of a point
  // input:  given point and number of nabors that are going to find
  // output: vertex id array and distance square array
  void nearestNabors(double x, double y, double z, 
    int*& vid, double*& dis2, double radius, int& k) {
    if (k>NbVertex()) k = NbVertex();
    if (pVertTree3==NULL) pVertTree3 = buildAnnKdTree();
    ANNdistArray dists = new ANNdist[k];	
    ANNidxArray idx = new ANNidx[k];
    ANNpoint p = annAllocPt(3);
    p[0] = x; p[1] = y; p[2] = z;
    // return the actual number of points nearby
    int kfound = pVertTree3->annkFRSearch(p, radius*radius, k, idx, dists, 0);
    k = min(k,kfound);
    vid = new int[k];
    dis2 = new double[k];
    for (int i=0; i<k; ++i) {
      vid[i] = idx[i];
      dis2[i] = dists[i];
    }
    annDeallocPt(p); delete [] dists; delete [] idx;
  }

  int nearestPnt(float x, float y, float z) {
    ANNdistArray dists = new ANNdist[1];	
    ANNidxArray idx = new ANNidx[1];
    ANNpoint p = annAllocPt(3);
    p[0] = x; p[1] = y; p[2] = z;
    if (pVertTree3==NULL) pVertTree3 = buildAnnKdTree();
    pVertTree3->annkSearch(p, 1, idx, dists, 0);
    int i = idx[0];
    delete [] idx, dists;
    annDeallocPt(p); 
    return i;
  }

  ///////////////////////////////////////////////////////////////////////
  // local reference frame by inverse-distance weighted 
  // Tombari et al. 2010 ECCV
  // the third eigen vector is the normal
  void localRF(double x, double y, double z, double radius, int k, 
    double* rfc) {
    // first find nearest neighbors inside the given radius
    int* vid;
    double* dis2;
    int i, j;
    nearestNabors(x,y,z,vid,dis2,radius,k);
    if (k<5) {
      cout<<"warning: num of nearest nabors less than 5"<<endl;
      rfc[0] = 1; rfc[1] = 0; rfc[2] = 0; 
      rfc[3] = 0; rfc[4] = 1; rfc[5] = 0;
      rfc[6] = 0; rfc[7] = 0; rfc[8] = 1;
      return;
    }
    // build covariance matrix
    double *covM[3], xdis[k], ydis[k], zdis[k], xd, yd, zd, dis, temp, sum=0;
    // Initialize covariance matrix
    covM[0] = new double[3];
    covM[1] = new double[3];
    covM[2] = new double[3];
    for (i=0; i<3; ++i)
      for (j=0; j<3; ++j) covM[i][j] = 0;
    for (i=0; i<k; ++i) {
      xdis[i] = xd = x-v(vid[i]).x();
      ydis[i] = yd = y-v(vid[i]).y();
      zdis[i] = zd = z-v(vid[i]).z();
      dis = sqrt(xd*xd+yd*yd+zd*zd);
      covM[0][0] += dis*xd*xd;
      covM[1][1] += dis*yd*yd;
      covM[2][2] += dis*zd*zd;
      temp = dis*xd*yd; covM[0][1] += temp; covM[1][0] += temp;
      temp = dis*xd*zd; covM[0][2] += temp; covM[2][0] += temp;
      temp = dis*yd*zd; covM[1][2] += temp; covM[2][1] += temp;
      sum += dis;
    } 
    double sum1 = 1.0/sum;
    for (i=0; i<3; ++i)
      for (j=0; j<3; ++j) covM[i][j] *= sum1;
    double eval[3], *evect[3];
    evect[0] = new double[3];
    evect[1] = new double[3];
    evect[2] = new double[3];
    // Diagonalization (eval = eigenvalues, evect = eigenvector)
    // - Eigenvalues and eigenvectors are sorted in decreasing order
    // - Eigenvectors are already normalized
    int resJ = vtkMath::Jacobi(covM, eval, evect);
    int plusNormal=0, plusTangentDirection1=0;
    double dp;
    for(i=0; i<k; ++i) {
      dp = xdis[i]*evect[0][0]+ydis[i]*evect[1][0]+zdis[i]*evect[2][0];	
      if (dp>=0) plusTangentDirection1++;
      dp = xdis[i]*evect[0][2]+ydis[i]*evect[1][2]+zdis[i]*evect[2][2];
      if (dp>=0) plusNormal++;
    }
    // PATCH: directions might still be ambiguous if point density is 
    // the same on both emispheres 
    bool isNormalStillAmbiguous=false, isTgDirStillAmbiguous=false;
    if(abs(plusTangentDirection1-k+plusTangentDirection1)==0) 
      isTgDirStillAmbiguous=true;
    if(abs(plusNormal-k+plusNormal)==0) isNormalStillAmbiguous=true;
    vector< pair<double, int> > vij_sorted;
    if (isNormalStillAmbiguous || isTgDirStillAmbiguous) {
      pair<double,int> tempP;
      for(i=0; i<k; ++i) {	
        tempP.first = sqrt(xdis[i]*xdis[i]+ydis[i]*ydis[i]+zdis[i]*zdis[i]);
        tempP.second = i;
        vij_sorted.push_back(tempP);
      }
      sort(vij_sorted.begin(),vij_sorted.end());
    }
    if( !isTgDirStillAmbiguous ) {	
      if (plusTangentDirection1 < k-plusTangentDirection1) {
        evect[0][0] *= -1;
        evect[1][0] *= -1;
        evect[2][0] *= -1;
      }
    }
    else {
      plusTangentDirection1=0;
      int pointsToDis = 5; ///std::min(valid_nn_points*2/2+1, 11);
      int nPointsToDis_half = 3;
      int indexToDis = k/2;
      double dotProduct;
      for (i=-pointsToDis/2; i<=pointsToDis/2; ++i) {
        dp = xdis[vij_sorted[indexToDis-i].second]*evect[0][0] + 
             ydis[vij_sorted[indexToDis-i].second]*evect[1][0] + 
             zdis[vij_sorted[indexToDis-i].second]*evect[2][0];
        if (dp>0) plusTangentDirection1++;	
      }		
      if (plusTangentDirection1<nPointsToDis_half){
        evect[0][0] *= -1;
        evect[1][0] *= -1;
        evect[2][0] *= -1;
      }
    }
    if(!isNormalStillAmbiguous) {
      if (plusNormal<k-plusNormal) {
        evect[0][2] *= -1;
        evect[1][2] *= -1;
        evect[2][2] *= -1;
      }
    }
    else{
      plusNormal=0;
      int pointsToDis = 5; ///std::min(valid_nn_points*2/2+1, 11);
      int nPointsToDis_half = 3;
      int indexToDis = k/2;
      double dotProduct;
      for (i=-pointsToDis/2; i<=pointsToDis/2; ++i){
        dp = xdis[vij_sorted[indexToDis-i].second]*evect[0][2] + 
             ydis[vij_sorted[indexToDis-i].second]*evect[1][2] + 
             zdis[vij_sorted[indexToDis-i].second]*evect[2][2];
        if (dotProduct>0) plusNormal++;	
      }
      if (plusNormal<nPointsToDis_half) {
        evect[0][2] *= -1;
        evect[1][2] *= -1;
        evect[2][2] *= -1;
      }	
    }
    rfc[6] = evect[0][2];
    rfc[7] = evect[1][2];
    rfc[8] = evect[2][2];
    rfc[0] = evect[0][0];
    rfc[1] = evect[1][0];
    rfc[2] = evect[2][0];
    // The last reference direction "n3" is the cross product of the other two
    vtkMath::Cross(&rfc[6], &rfc[0], &rfc[3]);
    for (i=0; i<3; ++i) {
      delete [] evect[i]; evect[i] = 0;
      delete [] covM[i]; covM[i] = 0;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  // normal estimation using local reference frames
  // vec3d n: regular outer normal for correction reason only
  void rfNormal(double x, double y, double z, double radius, int k, Vec3d n, 
    double& nx, double& ny, double& nz) {
    double* rfc = new double[9];
    localRF(x, y, z, radius, k, rfc);
    Vec3d nt(rfc[6], rfc[7], rfc[8]);
    if (Dot(nt,n)>0) {
      nx = rfc[6]; ny = rfc[7]; nz = rfc[8];
    } else {
      nx = -rfc[6]; ny = -rfc[7]; nz = -rfc[8];
    }
    delete [] rfc;
  }

  void rfNormal(double radius, int k) {
    double x, y, z, nx, ny, nz;
    computeNormal();
    for (int i=0; i<NbVertex(); ++i) {
      x = vTable[i].x();
      y = vTable[i].y();
      z = vTable[i].z();
      rfNormal(x,y,z,radius,k,vTable[i].getNormal(),nx,ny,nz);
      Vec3d nv(nx,ny,nz);
      vTable[i].setNormal(nv);
    }
  }
  
  ///////////////////////////////////////////////////////////////////////
  // unambiguous directions
  void pntDir(double x, double y, double z, 
    double& nx, double& ny, double&nz) {
    int i, j, nf = NbFace();
    double xx, yy ,zz, ar;
    nx = ny = nz = 0;
    const double onethird = 1.0/3;
    for (i=0; i<nf; ++i) {
      xx=0; yy=0; zz=0;
      ar = f(i).getArea();
      for (j=0; j<3; ++j) {
        xx += v(f(i).v(j)).x();
        yy += v(f(i).v(j)).y();
        zz += v(f(i).v(j)).z();
      }
      xx *= onethird; yy *= onethird; zz *= onethird;
      nx += (x-xx)*ar;
      ny += (y-yy)*ar;
      nz += (z-zz)*ar;
    }
    double len = nx*nx+ny*ny+nz*nz;
    if (len!=0) {
      double len1 = 1.0/sqrt(len);
      nx *= len1; ny *= len1; nz *= len1;
    }
    else {
      nx = 1; ny = 0; nz = 0;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  // local version of unambiguous directions
  void pntDir(double x, double y, double z, 
    double& nx, double& ny, double&nz, double ratioSupport) {
    int i, j, nf = NbFace();
    double xx, yy ,zz, ar, d2, radius = bBox.radius()*ratioSupport;
    nx = ny = nz = 0;
    const double onethird = 1.0/3;
    for (i=0; i<nf; ++i) {
      xx=0; yy=0; zz=0;
      ar = f(i).getArea();
      for (j=0; j<3; ++j) {
        xx += v(f(i).v(j)).x();
        yy += v(f(i).v(j)).y();
        zz += v(f(i).v(j)).z();
      }
      xx *= onethird; yy *= onethird; zz *= onethird;
      d2 = (x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz);
      if (d2>radius*radius) continue;
      nx += (x-xx)*ar;
      ny += (y-yy)*ar;
      nz += (z-zz)*ar;
    }
    double len = nx*nx+ny*ny+nz*nz;
    if (len!=0) {
      double len1 = 1.0/sqrt(len);
      nx *= len1; ny *= len1; nz *= len1;
    }
    else {
      nx = 1; ny = 0; nz = 0;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  // normal estimation by averaging face normals
  void computeNormal() {
    int i,v1,v2,v3;
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    for (i=0; i<NbVertex(); ++i) v(i).clearNormal();
    for (i=0; i<NbFace(); ++i) {
      v1=fTable[i].v1(); v2=fTable[i].v2(); v3=fTable[i].v3();
      x1=vTable[v1].x(); y1=vTable[v1].y(); z1=vTable[v1].z();
      x2=vTable[v2].x(); y2=vTable[v2].y(); z2=vTable[v2].z();
      x3=vTable[v3].x(); y3=vTable[v3].y(); z3=vTable[v3].z();
      // v1->v2
      Vec3d v12(x2-x1,y2-y1,z2-z1);
      // v1->v3
      Vec3d v13(x3-x1,y3-y1,z3-z1);
      // cross product
      Vec3d vc = Cross(v12,v13);
      Normalize(vc);
      fTable[i].setNormal(vc);
      vTable[v1].addNumFaceAroundByOne(); vTable[v1].addNormal(vc);
      vTable[v2].addNumFaceAroundByOne(); vTable[v2].addNormal(vc);
      vTable[v3].addNumFaceAroundByOne(); vTable[v3].addNormal(vc);
    }
    for (i=0; i<NbVertex(); ++i) {
      Vec3d nv = v(i).getNormal();
      nv /= v(i).getNumFaceAround();
      Normalize(nv);
      vTable[i].setNormal(nv);
    }
  }

  /**
   * average edge length =  mesh resolution
   */
  float avgEdgeLength() {
    float ael = 0;
    int ecount = 0, i, j, v1, v2;
    for (i=0; i<NbFace(); ++i) 
      for (j=0; j<3; ++j) {
        v1 = f(i).v(j);
        v2 = f(i).v((j+1)%3);
        ael += sqrt(v(v1).dis2(v(v2)));
        ecount++;
      }
    return ael/ecount;
  }

  /**
   * Add random noise to the mesh: ratio * mesh resolution (avg edge length)
   */
  void addRandomNoise(float ratio) {
    float meshRes = avgEdgeLength();
    for (int i=0; i<NbVertex(); ++i) {
      float x = (rand()%1000)/1000.0*ratio*meshRes;
      float y = (rand()%1000)/1000.0*ratio*meshRes;
      float z = (rand()%1000)/1000.0*ratio*meshRes;
      v(i).set(x+v(i).x(), y+v(i).y(), z+v(i).z());
    }
    // recompute normals
    computeNormal();
  }

  double fullArea() {
    double ar = 0;
    for (int i=0; i<NbFace(); ++i) ar += f(i).getArea();
    return ar;
  }

  double computeFaceArea() {
    for (int i=0; i<NbFace(); ++i) {
      Vec3d v21 = v(f(i).v1()).getV()-v(f(i).v2()).getV();
      Vec3d v31 = v(f(i).v1()).getV()-v(f(i).v3()).getV();
      Vec3d signedArea = Cross(v21,v31)*0.5;
      f(i).setArea(Length(signedArea));
    }
  }

  void draw(float alpha) {
    glEnable(GL_COLOR_MATERIAL);
    Vec3d nv,cv,nf;
    Face face;
    Vertex v1,v2,v3;
    glBegin(GL_TRIANGLES);
    for(int i=0; i<NbFace(); ++i) {
      face = f(i);
      v1 = v(face.v1()); 
      v2 = v(face.v2()); 
      v3 = v(face.v3());
      cv = v1.getColor();
      nv = v1.getNormal();
      nf = face.getNormal();
      //glNormal3f(nf.x,nf.y,nf.z);
      glNormal3f(nv.x,nv.y,nv.z);
      glColor4f(cv.x,cv.y,cv.z,alpha);
      glVertex3f(v1.x(),v1.y(),v1.z());
      nv = v2.getNormal();
      cv = v2.getColor();
      glNormal3f(nv.x,nv.y,nv.z);
      glColor4f(cv.x,cv.y,cv.z,alpha);
      glVertex3f(v2.x(),v2.y(),v2.z());
      nv = v3.getNormal();
      cv = v3.getColor();
      glNormal3f(nv.x,nv.y,nv.z);
      glColor4f(cv.x,cv.y,cv.z,alpha);
      glVertex3f(v3.x(),v3.y(),v3.z());
    }
    glEnd();
  }
  void drawPerFace(float alpha) {
    glEnable(GL_COLOR_MATERIAL);
    Vec3d nv,cv;
    Face face;
    Vertex v1,v2,v3;
    glBegin(GL_TRIANGLES);
    for(int i=0; i<NbFace(); ++i) {
      face = f(i);
      cv = face.getColor();
      nv = face.getNormal();
      v1 = v(face.v1()); 
      v2 = v(face.v2()); 
      v3 = v(face.v3());
      glNormal3f(nv.x,nv.y,nv.z);
      glColor4f(cv.x,cv.y,cv.z,alpha);
      glVertex3f(v1.x(),v1.y(),v1.z());
      glVertex3f(v2.x(),v2.y(),v2.z());
      glVertex3f(v3.x(),v3.y(),v3.z());
    }
    glEnd();
  }
  void drawComplexs() {
    glEnable(GL_COLOR_MATERIAL);
    glColor4f(0.3f,0.3f,0.3f,1);
    int v0, v1;
    float r = bBox.smax*0.01f;
    for (int i=0; i<cp->cellNum(); ++i) {
      cell* c = cp->get_cell(i);
      // draw 0d cells
      if (c->dim()==0) {
        int j = c->getVerInd(0);
        glPushMatrix();
        glTranslatef(v(j).x(),v(j).y(),v(j).z());
        //glutSolidSphere(r,16,16);
        glPopMatrix();
      }
      // draw 1d cells
      glDisable(GL_LIGHTING);
      glColor4f(1.0f,1.0f,1.0f,1);
      glLineWidth(5.0f);
      glBegin(GL_LINES);
      if (c->dim()==1) 
        for (int k=0; k<c->getEdgNum(); ++k) {
          ePair e = c->getEdg(k);
          v0 = e.v1; v1 = e.v2;
          glVertex3f(v(v0).x(),v(v0).y(),v(v0).z());
          glVertex3f(v(v1).x(),v(v1).y(),v(v1).z());
        }
      glEnd();
      glEnable(GL_LIGHTING);
      // draw 2d cells: center
      if (c->dim()==2) {
        int id = c->vid;
        glPushMatrix();
        glTranslatef(v(id).x(),v(id).y(),v(id).z());
        glutSolidSphere(r,16,16);
        glPopMatrix();
      }
    }
  }

  /**
   * For each vertex, build its nabor list, both face nabors and vertex
   * nabors.
   * Make sure all faces in this mesh have been oriented consistently before
   * calling this method.
   */
  void buildAdjacency() {
    int i, j, v1, v2, v3;
    for (i=0; i<NbVertex(); ++i) {
      v(i).clearFaceNabors();
      v(i).clearVertexNabors();
    }
    // add edge info
    eTable.clear();
    for (i=0; i<NbFace(); ++i) {
      for (j=0; j<3; ++j) {
        v1 = f(i).v(j);
        v2 = f(i).v((j+1)%3);
        // check redundancy
        int eInd=-1;
        Edge edg(v1,v2);
        for (int ii=0; ii<NbEdge(); ++ii) 
          if (eTable[ii].isSame(edg)) {
            eInd = ii;
            break;    
          }
        if (eInd==-1) {
          eTable.push_back(edg);
          eInd = eTable.size()-1;
        }
        v(v1).addEdgeNabor(eInd);
        v(v2).addEdgeNabor(eInd);
      }
    }
    // add face nabors for vertices
    for (i=0; i<NbFace(); ++i) 
      for (j=0; j<3; ++j) v(f(i).v(j)).addFaceNabor(i);
    // add vertex nabors for vertices
    for (i=0; i<NbVertex(); ++i) 
      for (j=0; j<v(i).NbFaceNabors(); ++j) {
        v1 = f(v(i).fNabor(j)).v1();
        v2 = f(v(i).fNabor(j)).v2();
        v3 = f(v(i).fNabor(j)).v3();
        // add a vertex neighbor to v(i) if 
        // the vertex is not as same as v(i) and 
        // the vertex has not been added to the neighbor vertex list before
        if (v1!=i && !v(i).hasVertexNabor(v1)) v(i).addVertexNabor(v1);
        if (v2!=i && !v(i).hasVertexNabor(v2)) v(i).addVertexNabor(v2);
        if (v3!=i && !v(i).hasVertexNabor(v3)) v(i).addVertexNabor(v3);
      }
    // add face nabors for faces
    for (i=0; i<NbFace(); ++i) 
      for (j=0; j<3; ++j) {
        // current vertex: f(i).v(j)
        // indices of other 2 vertices in this face are
        v1 = f(i).v((j+1)%3);        
        v2 = f(i).v((j+2)%3); 
        // loop over v1's nabor faces, find the one contains both v1 and v2,
        // but is NOT the current face
        for (int k=0; k<v(v1).NbFaceNabors(); ++k) 
          if (f(v(v1).fNabor(k)).contains(v2) && v(v1).fNabor(k)!=i) {
            f(i).setNabor(j,v(v1).fNabor(k));
            break;
          }
      }
  }

  /**
   * Find two face indices that share the input edge.
   * input:  the start and end point indices: a and b.
   * output: two face indices: f1 and f2.
   */
  void findSharedFaces(int a, int b, int& f1, int& f2) {
    int count = 0;
    // loop over a's face nabors
    for (int i=0; i<v(a).NbFaceNabors(); ++i) 
      // find out the ones that contain both a and b
      if (f(v(a).fNabor(i)).contains(b)) {
        if (count==0) {
          count++;
          f1 = v(a).fNabor(i);
        }
        else f2 = v(a).fNabor(i);
      }
  }

  /**
   * Find the third point in face f, when the first two points v1 and v2 
   * are given.
   **/
  int findThirdPoint(int v1, int v2, int f1) {
    for (int i=0; i<3; ++i) 
      if (f(f1).v(i)!=v1 && f(f1).v(i)!=v2) return f(f1).v(i);
  }

  /**
   * Generate monte carlo samples on the surface.
   * Input the expected number of samples, return the sample Point3 array.
   */
  vector<Vec3d> genRandomSamples(int n) {
    double s,t;
    vector<Vec3d> samples;
    if (n<1) return samples;
    computeFaceArea();
    // average area occupied by one sample
    double avgAr = fullArea()/n;
    // loop over faces to generate point samples
    for (int i=0; i<NbFace(); ++i) {
      Face fa = f(i);
      // three points
      Point3 p1(v(fa.v1()).x(), v(fa.v1()).y(), v(fa.v1()).z());
      Point3 p2(v(fa.v2()).x(), v(fa.v2()).y(), v(fa.v2()).z());
      Point3 p3(v(fa.v3()).x(), v(fa.v3()).y(), v(fa.v3()).z());
      // how many samples should be picked inside this triangle
      int nv = fa.getArea()/avgAr;
      // the probability of picking the last point
      double prob = fa.getArea()/avgAr-nv;
      // origin
      Vec3d o = v(fa.v1()).getV();
      // two edge coordinate axes
      Vec3d v12 = v(fa.v2()).getV()-o;
      Vec3d v13 = v(fa.v3()).getV()-o;
      Vec3d vs;
      // pick samples
      for (int j=0; j<nv; ++j) {
        // s\in(0,1)
        s = (rand()%1000)/1000.0;
        // t\in(0,1)
        t = (rand()%1000)/1000.0;
        if (s+t>1) {
          s = 1-s;
          t = 1-t;
        }
        vs = o+v12*s+v13*t;
        samples.push_back(vs);
      }
      // decide if I need to pick one more point
      double ss = (rand()%1000)/1000.0;
      if (ss<=prob) {
        s = (rand()%1000)/1000.0;
        t = (rand()%1000)/1000.0;
        if (s+t>1) {
          s = 1-s;
          t = 1-t;
        }
        vs = o+v12*s+v13*t;
        samples.push_back(vs);
      }
    }
    return samples;
  }

  /**
   * Generate monte carlo samples on the surface.
   * Input the expected number of samples, return the sample Point3 array.
   */
  void genRandomSamples(int n, vector<Point3>* psamp) {
    double s,t;
    psamp->clear();
    if (n<1) return;
    computeFaceArea();
    // average area occupied by one sample
    double avgAr = fullArea()/n;
    // loop over faces to generate point samples
    for (int i=0; i<NbFace(); ++i) {
      Face fa = f(i);
      // three points
      Point3 p1(v(fa.v1()).x(), v(fa.v1()).y(), v(fa.v1()).z());
      Point3 p2(v(fa.v2()).x(), v(fa.v2()).y(), v(fa.v2()).z());
      Point3 p3(v(fa.v3()).x(), v(fa.v3()).y(), v(fa.v3()).z());
      // how many samples should be picked inside this triangle
      int nv = fa.getArea()/avgAr+1;
      // the probability of picking the last point
      double prob = fa.getArea()/avgAr-nv;
      // origin
      Vec3d o = v(fa.v1()).getV();
      // two edge coordinate axes
      Vec3d v12 = v(fa.v2()).getV()-o;
      Vec3d v13 = v(fa.v3()).getV()-o;
      Vec3d vs;
      // pick samples
      for (int j=0; j<nv; ++j) {
        // s\in(0,1)
        s = (rand()%1000)/1000.0;
        // t\in(0,1)
        t = (rand()%1000)/1000.0;
        if (s+t>1) {
          s = 1-s;
          t = 1-t;
        }
        vs = o+v12*s+v13*t;
        // estimate normal: linear interpolation use barycentric coordinates
        double l1, l2, l3;
        Point3 p(vs.x,vs.y,vs.z);
        baryCoord(p,p1,p2,p3,l1,l2,l3);
        Vec3d no = v(fa.v1()).getN()*l1+v(fa.v2()).getN()*l2+
          v(fa.v3()).getN()*l3;
        //p.set(vs.x,vs.y,vs.z,no.x,no.y,no.z);
        p.set(vs.x,vs.y,vs.z,fa.nx(),fa.ny(),fa.nz());
        psamp->push_back(p);
      }
      // decide if I need to pick one more point
      double ss = (rand()%1000)/1000.0;
      if (ss<=prob) {
        s = (rand()%1000)/1000.0;
        t = (rand()%1000)/1000.0;
        if (s+t>1) {
          s = 1-s;
          t = 1-t;
        }
        vs = o+v12*s+v13*t;
        // estimate normal: linear interpolation use barycentric coordinates
        double l1, l2, l3;
        Point3 p(vs.x,vs.y,vs.z);
        baryCoord(p,p1,p2,p3,l1,l2,l3);
        Vec3d no = v(fa.v1()).getN()*l1+v(fa.v2()).getN()*l2+
          v(fa.v3()).getN()*l3;
        p.set(vs.x,vs.y,vs.z,fa.nx(),fa.ny(),fa.nz());
        //p.set(vs.x,vs.y,vs.z,no.x,no.y,no.z);
      }
    }
  }

  /************************************************************************/
  // consistently oriente faces
  void push(int i) {
    if (flags[i]) return;
    assert(sp<NbFace());
    stack[sp++] = i;
  }

  int pop() {
    assert(sp>0);
    return stack[--sp];
  }

  void reorient(int i) {
    reorients++;
    f(i).reorient();
  }

  int which(int i, int t) {
    if (i==f(t).v1()) return 0;
    if (i==f(t).v2()) return 1;
    assert(i==t.v3());
    return 2;
  }

  int isalong(int v1, int v2, int t) {
    return which(v2,t)==(which(v1,t)+1)%3;
  }

  void _traverse() {
    int i,j,k;
    while(sp>0) {
      i = pop();
      flags[i] = 1;
      for (j=0; j<3; ++j) {
        // find neighbor face which opposite vertex j
        k = f(i).getNabor(j);
        assert(k!=i);
        if (flags[k]) continue;
        /* fix the orientation of k if necessary */
        if (isalong(f(i).v((j+1)%3),f(i).v((j+2)%3),k)) reorient(k);
        push(k);
      }
    }
  }

  void traverse(int i) {
    assert(!sp);
    push(i);
    _traverse();
  }

  void orientFaces() {
    int size = NbFace();
    stack = (int*)malloc(size*sizeof(int));;
    sp = 0;
    reorients = 0;
    flags = (unsigned char*)malloc(size*sizeof(unsigned char));
    memset(flags,0,size*sizeof(unsigned char));
    for (int i=0; i<size; ++i)
      if (!flags[i]) {
        traverse(i);
      }
    free(flags);
    free(stack);
    cout<<reorients<<" orient faces completed"<<endl;
  }

};
}

#endif
