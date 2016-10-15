
#ifndef __MESH_BASE_H
#define __MESH_BASE_H

#include <global.h>
#include <vector>

/* ------------------------------------------------------ */

class mesh_element {

 public:

  int ID;
  int dimension;

  int faces;
  mesh_element** face;   // ALL faces, in CCW order for 2D face 
  // for a 2D face 0D (vertices) are odd, 1D (edges) are even

  int cofaces;
  mesh_element** coface; // ALL cofaces; for a vertex, in order around it

  mesh_element ( int id, int dim );
  ~mesh_element();

  void set_faces ( int n, mesh_element **f );
  void set_cofaces ( int n, mesh_element **cf );
  bool isboundary();
  int find_face_index ( mesh_element *e );  // which face is this?

  bool isisolatedvertex();

  void print_out();
};

/* ------------------------------------------------------ */

// assumes manifold mesh, with planar 2D cells

class mesh_base {

  std::vector<std::vector<int>*> *workspace;

 protected:

  int vtcs, es, fcs;
  std::vector<mesh_element*> mel;  // mesh elements in decreasing dimension order (2,1,0)

 public:
  void add_2Dmel ( std::vector<int> *verts );   // adds a 2D mesh element
  void finalize();    // finalizes the datastructure; 
                      // in particular, computes faces, cofaces etc
  mesh_base ( );      
  ~mesh_base();

  mesh_element * get ( int i ); // get mesh element i
  mesh_element * getedge ( int i );
  mesh_element * getvertex ( int i );
  mesh_element * getface ( int i ); 

  void print_out();

  int vertices();
  int edges();
  int faces();
  int mesh_elements();
};

/* ------------------------------------------------------ */

#endif
