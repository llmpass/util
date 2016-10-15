
#include <vector>
#include <ostream>
#include <mesh.h>

using namespace std;

#pragma once

/* ------------------------------------------ */

class complexs;
class cell;
class ePair {
  public:
  ePair(int i1, int i2) : v1(i1), v2(i2) {}
  int v1, v2;
};

class cell {
  
  friend class ePair;
  friend class complexs;
  friend ostream & operator<< ( std::ostream &o, const cell &c );
  
  // counters used to generate IDs of faces of different dimensions
  static int curid0;
  static int curid1;
  static int curid2;
  
  int ID;
  bool isalive;               // is in complex?
  int d;                      // dimension

  vector<int> coef;      // coefficients of the boundary operator
  vector<cell*> f;       // faces : \sum coef[i]*f[i] is the boundary
  vector<cell*> c;       // cofaces: 0d cells have 1d cofaces and 1d cells
                         // have 2d cofaces
  vector<int> triInd;    // indices of contained triangles
  vector<ePair> edg;    // contained edges
  vector<int> verInd;    // indices of contained vertices

  // remove cell from (co)face list
  // returns items removed
  int remove_face ( const cell *const o );
  int remove_coface ( const cell *const o );

  // consistency check stuff...
  bool has_same_cofaces ( const cell *const o ) const;

  // how many cofaces equal to o
  int has_coface ( const cell *const o ) const;   
  int has_face ( const cell *const o ) const;

  // same, but algebraic
  int a_has_face ( const cell *const o ) const;

  // de-associate the cell
  void unface();
  void uncoface();

  double siz;

 public:
  // ancestor vertex id
  int vid;
  // average geodesic distance
  float agd;
  // number of faces
  int faces() const;
  // number of cofaces
  int cofaces() const;
  // dimension (0,1 or 2)
  int dim() const;
  operator bool() const;  // is the cell alive?

  // return face, coface or coef #i
  cell * const &face ( int i ) const;
  int coefficient ( int i ) const;
  cell * const &coface ( int i ) const;
  cell *&face ( int i );
  cell *&coface ( int i );
  vector<int>& getTriInd() {return triInd;}
  vector<ePair>& getEdg() {return edg;}
  void addTriInd(int triId) {triInd.push_back(triId);}
  void addVerInd(int verId) {verInd.push_back(verId);}
  void addEdg(ePair& e) {edg.push_back(e);}
  int getTriInd(int i) {return triInd[i];} 
  int getVerInd(int i) {return verInd[i];} 
  ePair& getEdg(int i) {return edg[i];}
  int getEdgNum() {return edg.size();}
  int getTriNum() {return triInd.size();}
  vector<int> getNabor2dCellId();

  double size() const;
  int id() const;

  // these methods also update nf's data!
  void add_face ( cell *nf, int cf );

  // merge cell c1 into this one
  //   c1 and c2 have to be adjacent
  //   c1 dies, the calling one will still be there, with updated boundary
  bool merge ( cell *c1 );

  cell ( int dd, double sz );

  // invariant/consistency check; returns false on fail
  bool check() const;

  // resets counters; always call before starting to construct 
  // a new complex!
  static void reset_counters();
};

/* ------------------------------------------ */

class complexs {

  friend ostream & operator<< ( std::ostream &o, const complexs &c );

  vector<cell*> c;

 public:
  
  // note: don't start constructing a new complex before finishing other:
  // do it one at a time

  int cells (int dim, int*& map1, int*& map2) const;

  // how many cells, dead or alive
  int cells_dead_or_alive ( ) const;

  // get cell #i (note it can be a dead one too)
  cell *get_cell ( int i ) const;

  // add new cell to complex
  int add_cell ( cell * nc );

  // make c1 face of c0, pointer and integer index versions
  void makeface ( cell *c0, cell *c1, int cf );
  void makeface ( int c0, int c1, int cf );

  complexs();
  complexs( const char *name );   // read complex from a .t file
  complexs( mesh *m );           // build a complex from a mesh

  ~complexs();

  // invariant/consistency check; returns false on fail
  bool check();

  // remove dead cells; note that this reindexes/changes IDs of all
  // alive cells...
  void remove_dead();
  int cellNum() {return c.size();}
};

/* ------------------------------------------ */

ostream & operator<< ( ostream &o, const cell &c );
ostream & operator<< ( ostream &o, const complexs &c );

/* ------------------------------------------ */
