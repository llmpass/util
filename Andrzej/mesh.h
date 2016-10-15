
#ifndef __MESH_H
#define __MESH_H

#include <global.h>
#include <fstream>
#include <mesh_base.h>
#include <vec3d.h>

#include <GLUT/glut.h>

/* ------------------------------------------------------ */

class mesh : public mesh_base {

  void read_from_file ( std::ifstream &ifs, char format );
  // format: 't' - .t file
  // format: 'n' - polygonal mesh file
  // like .t but faces terminated with -1
  // file in n format starts with 'n'


  vec3dd evec ( const mesh_element *f, int i, int j );

 protected:

  std::ifstream ifs;

 public:

  vec3dd *v;  // vertex coordinates
  vec3dd *n;  // unit normals for the faces
  int* group; // group index of the faces

  mesh() {v=n=NULL;} // blank mesh
  mesh ( const char *name ); 
  ~mesh();

  void print_out();
  void compute_normals();
  // returns point on edge eid 
  vec3dd edgepoint ( int eid, double t );  
  vec3dd vertex ( int i ); // vertex coordinate
  vec3dd normal ( int i ); // face normal
  void draw(vec3dd c);
};

/* ------------------------------------------------------ */

#endif
