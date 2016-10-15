#include <mesh.h>
#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

/* ------------------------------------------------------ */

mesh::mesh ( const char *name ) : mesh_base(), ifs(name)
{
  if (!ifs) 
    {
      cout << "Can't open " << name << endl;
      assert(ifs);
    }

  // detect format
  char c = ifs.get();
  char format;

  if (c=='n')
    format = 'n';
  else
    {
      ifs.putback(c);
      format = 't';
    }

  read_from_file(ifs,format);

  compute_normals();
}

/* ------------------------------------------------------ */

void mesh::read_from_file ( ifstream &ifs, char format )
{
  switch(format)
    {
    case 't':
      {
	int vs,ts,i;
	ifs >> ts >> vs;
	for ( i=0; i<ts; i++ )
	  {
	    int a,b,c;
	    ifs >> a >> b >> c;
	    vector<int> *t = new vector<int>;
	    t->push_back(a);
	    t->push_back(b);
	    t->push_back(c);
	    add_2Dmel(t);
	  }
	finalize();

	if (vs!=vtcs)
	  {
	    cout << "Warning: nused vertices detected" << endl;
	  }

	v = new vec3dd[vtcs];
	for ( i=0; i<vtcs; i++ )
	  ifs >> v[i][0] >> v[i][1] >> v[i][2];
      }
      break;

    case 'n':
      {
	int vs,fs,i;
	ifs >> fs >> vs;
	for ( i=0; i<fs; i++ )
	  {
	    int a;
	    vector<int> *t = new vector<int>;
	    while(1)
	      {
		ifs >> a;
		if (a==-1)
		  break;
		t->push_back(a);
	      }
	    add_2Dmel(t);
	  }
	finalize();

	if (vs!=vtcs)
	  {
	    cout << "Warning: nused vertices detected" << endl;
	  }

	v = new vec3dd[vtcs];
	for ( i=0; i<vtcs; i++ )
	  ifs >> v[i][0] >> v[i][1] >> v[i][2];
      }
    }
}

/* ------------------------------------------------------ */

vec3dd mesh::evec ( const mesh_element *f, int i, int j )
{
  return v[f->face[j]->ID]-v[f->face[i]->ID];
}

/* ------------------------------------------------------ */

vec3dd mesh::edgepoint ( int eid, double t )
{
  mesh_element *e = getedge(eid);
  return (1-t)*v[e->face[0]->ID]+t*v[e->face[1]->ID];
}

/* ------------------------------------------------------ */

vec3dd mesh::vertex ( int i )
{
  return v[i];
}

/* ------------------------------------------------------ */

vec3dd mesh::normal ( int i )
{
  return n[i];
}

/* ------------------------------------------------------ */

void mesh::compute_normals()
{
  n = new vec3dd[fcs];
  for ( int i=0; i<fcs; i++ )
    {
      vec3dd n0(0,0,0);
      const mesh_element *f = get(i);
      for ( int j=1; j<f->faces-4; j+=2 )
	n0 += evec(f,j,j+2)^evec(f,j,j+4);
      n0.normalize();
      n[i] = n0;
    }
}

/* ------------------------------------------------------ */

void mesh::print_out()
{
  int i;
  mesh_base::print_out();
  cout << " ---------------------- " << endl;
  cout << "Vertex coordinates" << endl;
  for ( i=0; i<vtcs; i++ )
    cout << i << " " << v[i] << endl;
  cout << " ---------------------- " << endl;
  cout << "Normals" << endl;
  for ( i=0; i<fcs; i++ )
    cout << i << " " << n[i] << endl;
}

/* ------------------------------------------------------ */

mesh::~mesh()
{
  if (v) delete[] v;
  if (n) delete[] n;
  v = n = NULL;
}

/* ------------------------------------------------------ */

// something added by Luming Liang

void mesh::draw(vec3dd c) {
  glEnable(GL_COLOR_MATERIAL);
  int i, j, v1, v2, v3;
  glColor3f(c[0],c[1],c[2]);
  glBegin(GL_TRIANGLES);
  for (i=0; i<fcs; ++i) {
    const mesh_element *f = getface(i);
    vec3dd nor = normal(i);
    glNormal3f(nor[0],nor[1],nor[2]);
    // the vertices are indexed with odd numbers
    int v1 = f->face[1]->ID;
    int v2 = f->face[3]->ID;
    int v3 = f->face[5]->ID;
    vec3dd vv1 = vertex(v1);
    glVertex3f(vv1[0],vv1[1],vv1[2]);
    vec3dd vv2 = vertex(v2);
    glVertex3f(vv2[0],vv2[1],vv2[2]);
    vec3dd vv3 = vertex(v3);
    glVertex3f(vv3[0],vv3[1],vv3[2]);
  }
  glEnd();
}
