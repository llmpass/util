
#include <cassert>
#include <complexs.h>
#include <mesh.h>

using namespace std;

/* ------------------------------------------ */

int cell::curid0 = 0;
int cell::curid1 = 0;
int cell::curid2 = 0;

/* ------------------------------------------ */

void cell::reset_counters()
{
  curid0 = curid1 = curid2 = 0;
}

/* ------------------------------------------ */

cell::cell ( int dd, double sz )  : siz(sz), d(dd), f(), c(), isalive(true) 
{
  switch(dd)
    {
    case 0:
      ID = curid0++;
      break;
    case 1:
      ID = curid1++;
      break;
    case 2:
      ID = curid2++;
      break;
    default:
      assert(0);
    }
}

/* ------------------------------------------ */

void cell::add_face ( cell *nf, int cf )
{
  assert(nf->d==d-1);
  f.push_back(nf);
  coef.push_back(cf);
  assert(f.size()==coef.size());
  nf->c.push_back(this);
}

/* ------------------------------------------ */

/*
void cell::add_coface ( cell *nf  )
{
  assert(nf->d==d+1);
  c.push_back(nf);
  nf->f.push_back(this);
}
*/

/* ------------------------------------------ */

double cell::size() const
{
  return siz;
}

/* ------------------------------------------ */

int cell::faces() const
{
  assert(f.size()==coef.size());
  return f.size();
}

/* ------------------------------------------ */

int cell::cofaces() const
{
  return c.size();
}

/* ------------------------------------------ */

cell * &cell::face ( int i ) 
{
  assert(i>=0 && i<faces());
  return f[i];
}

/* ------------------------------------------ */

int cell::coefficient ( int i ) const
{
  assert(i>=0 && i<faces());
  return coef[i];
}

/* ------------------------------------------ */

cell * &cell::coface ( int i ) 
{
  assert(i>=0 && i<cofaces());
  return c[i];
}

/* ------------------------------------------ */

cell * const &cell::face ( int i ) const
{
  assert(i>=0 && i<faces());
  return f[i];
}

/* ------------------------------------------ */

cell * const &cell::coface ( int i ) const 
{
  assert(i>=0 && i<cofaces());
  return c[i];
}

/* ------------------------------------------ */

bool cell::merge ( cell *c1 ) {
  if (d==c1->d && d==1) {
      // 1D edge merger

      cell *sharedv;  // the shared vertex
      assert(faces()==2 && c1->faces()==2);

      if (face(0)==c1->face(0) || face(0)==c1->face(1))
	sharedv = face(0);
      if (face(1)==c1->face(0) || face(1)==c1->face(1))
	sharedv = face(1);

      // is the merger legal? only if the other vertices are 
      // not the same

      cell *otherv0 = face(face(0)==sharedv ? 1 : 0);
      cell *otherv1 = c1->face(c1->face(0)==sharedv ? 1 : 0);
      assert(otherv0!=sharedv);
      assert(otherv1!=sharedv);
      if (otherv0==otherv1)
	return false; // can't merge

      // ready to merge now !

      // this should be true for closed meshes
      assert(has_same_cofaces(c1));

      // cofaces are fine -- fix the face list
      // shared vertex needs to be replaced
      if (sharedv==face(0))
	face(0) = otherv1;
      else
	{
	  assert(sharedv==face(1));
	  face(1) = otherv1;
	}

      siz += c1->siz;
      // add agd of c1 into the current 1d cell
      agd = (c1->siz*c1->agd+siz*agd)/(c1->siz+siz);

      // finally, c1 should die now.
      c1->isalive = false;
      // remove all associacions with other faces
      c1->unface();
      c1->uncoface();
      otherv1->c.push_back(this);

      // and so does the vertex...
      sharedv->isalive = false;
      sharedv->uncoface();
      sharedv->unface();

      // merge edge indices contained in c1 to the current cell
      for (int ii=0; ii<c1->getEdgNum(); ++ii) 
        edg.push_back(c1->getEdg(ii));
      return true;
    }

  if (d==c1->d && d==2)
    {
      // 2D merger
      
      if (this==c1) return false;  // can't merge with itself
      assert(!cofaces());
      assert(!c1->cofaces());

      vector<cell*> newf;
      vector<cell*> newf2;
      vector<cell*> toremove;
      vector<int> newcf;

      int flip = -1;

      for ( int ii=0; ii<f.size(); ii++ )
	{
	  cell *i = f[ii];
	  int mme = has_face(i);
	  int mc1 = c1->has_face(i);
	  int amme = a_has_face(i);
	  int amc1 = c1->a_has_face(i);
	  assert(mme>=0 && mme<=1);
	  assert(mc1>=0 && mc1<=1);
	  assert(amme>=-1 && amme<=1);
	  assert(amc1>=-1 && amc1<=1);
	  assert((amme+2)%2==mme);
	  assert((amc1+2)%2==mc1);
	  if ((mme+mc1) % 2)
	    {
	      newf.push_back(i);
	      newcf.push_back(coef[ii]);
	    }
	  else
	    {
	      toremove.push_back(i);
	      if (amme+amc1)
		flip = 1;
	      else
		flip = 0;
	    }
	  int test1 = i->remove_coface(c1);
	  int test2 = c1->remove_face(i);
	  assert(test1<=1 && test1>=0);
	  assert(test2<=1 && test2>=0);
	  assert(test1==test2);
	}

      assert(flip>=0);

      for ( int ii=0; ii<c1->f.size(); ii++ )
	{
	  cell *i = c1->f[ii];
	  assert(!has_face(i));
	  int mc1 = c1->has_face(i);
	  assert(mc1==1);	  
	  newf.push_back(i);
	  newcf.push_back(flip ? -c1->coef[ii] : c1->coef[ii]);
	  newf2.push_back(i);
	}  
      // merge triangle indices contained in c1 to the current cell
      vector<int> triInd1 = c1->getTriInd();
      for (int ii=0; ii<triInd1.size(); ++ii) 
        triInd.push_back(triInd1[ii]);
      // update average AGD of the current cell
      agd = (c1->siz*c1->agd+siz*agd)/(c1->siz+siz);

      f = newf;
      coef = newcf;

      siz += c1->siz;

      // c1 should die now.
      c1->isalive = false;
      // remove all associacions with other faces
      c1->unface();
      c1->uncoface();

      for ( int ii=0; ii<toremove.size(); ii++ )
	{
	  cell *i = toremove[ii];
	  i->isalive = false;
	  i->unface();
	  i->uncoface();
	}

      for ( int ii=0; ii<newf2.size(); ii++ )
	{
	  cell *i = newf2[ii];
	  i->c.push_back(this);
	}

      // and finally, merge edges on the boundary of the new face...
      for ( int i=0; i<faces(); )
	{
	  cell *f0 = face(i);
	  assert(f0->faces()==2);
	  assert(f0->d==1);

	  bool mgd = false;

	  for ( int j=0; j<2; )
	    {
	      cell *c = f0->face(j);
	      assert(c->cofaces()>=2);
	      if (c->cofaces()!=2)
		{
		  ++j;
		  continue;
		}

	      cell *f1 = 0;
	      if (c->coface(0)==f0)
		f1 = c->coface(1);
	      else
		{
		  assert(c->coface(1)==f0);
		  f1 = c->coface(0);
		}
	      assert(f1);

	      bool isok = f0->merge(f1);
	      mgd = mgd | isok;

	      if (!isok) ++j;
	    }

	  if (!mgd) i++;
          else i=0;
	}

      return true;
    }

  return false;
}

/* ------------------------------------------ */

bool cell::has_same_cofaces ( const cell *const o ) const 
{
  for ( int ii=0; ii<c.size(); ii++ )
    {
      cell *i = c[ii];
      if (has_coface(i)!=o->has_coface(i))
	return false;
    }
  for ( int ii=0; ii<c.size(); ii++ )
    {
      cell *i = c[ii];
      if (has_coface(i)!=o->has_coface(i))
	return false;
    }
  return true;
}

/**
 * Return cell Ids of naboring 2d cells of a 0d cell.
 */
vector<int> cell::getNabor2dCellId() {
  vector<int> c2dAroundId;
  int i,j,k;
  // find 2d cells that locate around this 0d cell
  // the cofaces of a 0d cell is 1d, the cofaces of a 1d cell is 2d
  // loop over 1d cofaces of the current 0d cell
  for (i=0; i<cofaces(); ++i) {
    cell* c1d = coface(i);
    // store the unstored 2d cofaces' Id of the 1d coface cf1[i]
    for (k=0; k<c1d->cofaces(); ++k) {
      cell* c2d = c1d->coface(k);
      bool hasStored = false;
      for (j=0; j<c2dAroundId.size(); ++j) 
        if (c2dAroundId[j]==c2d->id()) {
          hasStored = true;
          break;
        }
      if (!hasStored) c2dAroundId.push_back(c2d->id());
    }
  }
  return c2dAroundId;
}

/* ------------------------------------------ */

int cell::has_coface ( const cell *const o ) const
{
  int res = 0;
  for ( int ii=0; ii<c.size(); ii++ )
    {
      cell *i = c[ii];
      if (i==o) ++res;
    }
  return res;
}

/* ------------------------------------------ */

int cell::has_face ( const cell *const o ) const
{
  int res = 0;
  for ( int ii=0; ii<f.size(); ii++ )
    {
      cell *i = f[ii];
      if (i==o) ++res;
    }
  return res;
}

/* ------------------------------------------ */

int cell::a_has_face ( const cell *const o ) const
{
  int res = 0;
  for ( int i=0; i<faces(); i++ )
    if (f[i]==o) res+=coef[i];
  return res;
}

/* ------------------------------------------ */

int cell::remove_face ( const cell *const o )
{
  int res = 0;
  for ( int i=0; i<faces(); )
    if (f[i]==o)
      {
	f[i] = f[faces()-1];
	coef[i] = coef[faces()-1];
	f.pop_back();
	coef.pop_back();
	res++;
      }
    else
      i++;
  return res;
}

/* ------------------------------------------ */

int cell::remove_coface ( const cell *const o )
{
  int res = 0;
  for ( int i=0; i<cofaces(); )
    if (c[i]==o)
      {
	c[i] = c[cofaces()-1];
	c.pop_back();
	res++;
      } 
    else
      i++;
  return res;
}

/* ------------------------------------------ */

void cell::unface() {
  for (int ii=0; ii<f.size(); ++ii) {
    cell *i = f[ii];
    i->remove_coface(this);
    if (i->d==0 && i->c.size()==0) i->isalive = false;
  }
  f.clear();
  coef.clear();
}

/* ------------------------------------------ */

void cell::uncoface()
{
  for ( int ii=0; ii<c.size(); ii++ )
    {
      cell *i = c[ii];
      i->remove_face(this);
    }
  c.clear();
}

/* ------------------------------------------ */

int cell::id() const
{
  return ID;
}

/* ------------------------------------------ */

bool cell::check() const
{
  if (!isalive)
    {
      return !(faces() || cofaces());
    }

  if (d==2)
    {
      if (cofaces()) return false;
      for ( int ii=0; ii<f.size(); ii++ )
	{
	  cell *i = f[ii];
	  if (has_face(i)!=1) return false;
	}

      // check the boundary operator....
      for ( int ii=0; ii<f.size(); ii++ )
	{
	  cell *i = f[ii];
	  for ( int jj=0; jj<i->f.size(); jj++ )
	    {
	      cell *j = i->f[jj];
	      int total = 0;
	      for ( int ii = 0; ii<f.size(); ii++ )
		for ( int jj = 0; jj<f[ii]->faces(); jj++ )
		  if (f[ii]->f[jj]==j)
		  total += coef[ii]*f[ii]->coef[jj];
	      if (total!=0)
		return false;
	    }
	}
    }

  if (d==1)
    {
      if (faces()!=2) return false;
      if (face(0)==face(1)) return false;
      if (cofaces()!=2) return false;
      if (coface(0)==coface(1)) return false;
    }

  if (d==0)
    {
      if (faces()) return false;
      for ( int ii=0; ii<c.size(); ii++ )
	{
	  cell *i = c[ii];
	  if (has_coface(i)!=1) return false;
	}
    }

  return true;
}

/* ------------------------------------------ */

ostream & operator<< ( ostream &o, const cell &c )
{
  o << "dimension: " << c.d << ", ID: " << c.ID << " ";
  if (!c.isalive) 
    o << "X" << endl;
  else
    o << endl;
  o << "    size: " << c.size() << endl;
  o << "    agd: " << c.agd << endl;
  o << "    faces: ";
  for ( int i=0; i<c.f.size(); i++ )
    o << "(" << c.f[i]->d << "," << c.f[i]->ID << ") ";
  o << endl;
  o << "    coefs: ";
  for ( int i=0; i<c.coef.size(); i++ )
    o << c.coef[i] << " ";
  o << endl;
  o << "    cofaces: ";
  for ( int i=0; i<c.c.size(); i++ )
    o << "(" << c.c[i]->d << "," << c.c[i]->ID << ") ";
  o << endl;
  return o;
}

/* ------------------------------------------ */
/* ------------------------------------------ */

complexs::complexs() : c() {
  cell::reset_counters();
}

/* ------------------------------------------ */

complexs::complexs ( mesh *m ) {
  cell::reset_counters();

  vector<int> f0;
  vector<int> f1;
  vector<int> f2;

  f0.reserve(m->vertices());
  f1.reserve(m->edges());
  f2.reserve(m->faces());

  for (int i=0; i<m->faces(); i++) {
    assert(m->getface(i)->ID==i);
    vec3dd area = 0;
    for (int j=1; j<m->getface(i)->faces; j+=2) {
      int jj  = (j+2)%m->getface(i)->faces;
      vec3dd a1 = m->vertex(m->getface(i)->face[1]->ID);
      vec3dd a2 = m->vertex(m->getface(i)->face[j]->ID);
      vec3dd a3 = m->vertex(m->getface(i)->face[jj]->ID);
      area += (a2-a1)^(a3-a1);
    }
    area *= 0.5;
    cell* c2 = new cell(2,area.norm());
    c2->addTriInd(i);
    c2->vid = m->getface(i)->face[1]->ID;
    f2[i] = add_cell(c2);
  }
  
  for (int i=0; i<m->edges(); i++) {
      assert(m->getedge(i)->ID==i);
      vec3dd v1 = m->vertex(m->getedge(i)->face[0]->ID)-
	m->vertex(m->getedge(i)->face[1]->ID);
      // find an edge that pointing from m->getedge(i)->face[0]->ID to 
      // m->getedge(i)->face[1]->ID
      cell* c1 = new cell(1,v1.norm());
      ePair e(m->getedge(i)->face[0]->ID,m->getedge(i)->face[1]->ID);
      c1->addEdg(e);
      c1->vid = m->getedge(i)->face[0]->ID;
      f1[i] = add_cell(c1);
    }

  for ( int i=0; i<m->vertices(); i++ )
    {
      assert(m->getvertex(i)->ID==i);
      cell* c0 = new cell(0,0.0);
      c0->addVerInd(i);
      c0->vid = i;
      f0[i] = add_cell(c0);
    }

  for ( int i=0; i<m->faces(); i++ )
    {
      mesh_element *me = m->getface(i);
      for ( int j=0; j<me->faces/2; j++ )
	makeface(f2[i],f1[me->face[2*j]->ID],
		 (me->face[2*j]->face[1]==me->face[2*j+1]) ? 1 : -1);
    }

  for ( int i=0; i<m->edges(); i++ )
    {
      mesh_element *me = m->getedge(i);
      for ( int j=0; j<me->faces; j++ )
	makeface(f1[i],f0[me->face[j]->ID],1-2*j);
    }
}

/* ------------------------------------------ */

complexs::complexs ( const char *name )
{
  cell::reset_counters();

  mesh m(name);

  vector<int> f0;
  vector<int> f1;
  vector<int> f2;

  f0.reserve(m.vertices());
  f1.reserve(m.edges());
  f2.reserve(m.faces());

  for ( int i=0; i<m.faces(); i++ )
    {
      assert(m.getface(i)->ID==i);
      vec3dd area = 0;
      for ( int j=1; j<m.getface(i)->faces; j+=2 )
	{
	  int jj  = (j+2)%m.getface(i)->faces;
	  vec3dd a1 = m.vertex(m.getface(i)->face[1]->ID);
	  vec3dd a2 = m.vertex(m.getface(i)->face[j]->ID);
	  vec3dd a3 = m.vertex(m.getface(i)->face[jj]->ID);
	  area += (a2-a1)^(a3-a1);
	}
      area *= 0.5;
      f2[i] = add_cell(new cell(2, area.norm()));
    }
  
  for ( int i=0; i<m.edges(); i++ )
    {
      assert(m.getedge(i)->ID==i);
      vec3dd v1 = m.vertex(m.getedge(i)->face[0]->ID)-
	m.vertex(m.getedge(i)->face[1]->ID);
      f1[i] = add_cell(new cell(1, v1.norm()));
    }

  for ( int i=0; i<m.vertices(); i++ )
    {
      assert(m.getvertex(i)->ID==i);
      f0[i] = add_cell(new cell(0,0.0));
    }

  for ( int i=0; i<m.faces(); i++ )
    {
      mesh_element *me = m.getface(i);
      for ( int j=0; j<me->faces/2; j++ )
	makeface(f2[i],f1[me->face[2*j]->ID],
		 (me->face[2*j]->face[1]==me->face[2*j+1]) ? 1 : -1);
    }

  for ( int i=0; i<m.edges(); i++ )
    {
      mesh_element *me = m.getedge(i);
      for ( int j=0; j<me->faces; j++ )
	makeface(f1[i],f0[me->face[j]->ID],1-2*j);
    }
}

/* ------------------------------------------ */

complexs::~complexs()
{
  for ( int i=0; i<c.size(); i++ )
    delete c[i];
  c.clear();
}

/* ------------------------------------------ */

int complexs::add_cell ( cell * nc )
{
  c.push_back(nc);
  return c.size()-1;
}

/* ------------------------------------------ */

void complexs::makeface ( cell *c0, cell *c1, int cf )
{
  c0->add_face(c1,cf);
}

/* ------------------------------------------ */

void complexs::makeface ( int c0, int c1, int cf )
{
  makeface(get_cell(c0),get_cell(c1),cf);
}

/* ------------------------------------------ */
/*
void complexs::makecoface ( cell *c0, cell *c1 )
{
  c0->add_coface(c1);
}
*/
/* ------------------------------------------ */
 /*
void complexs::makecoface ( int c0, int c1 )
{
  makecoface(get_cell(c0),get_cell(c1));
}
*/
/* ------------------------------------------ */

cell *complexs::get_cell ( int i ) const
{
  assert(i>=0 && i<c.size());
  return c[i];
}

/* ------------------------------------------ */

bool complexs::check()
{
  for ( int i=0; i<c.size(); i++ )
    if (!c[i]->check()) return false;
  return true;
}

/* ------------------------------------------ */

ostream & operator<< ( ostream &o, const complexs &c )
{
  for ( int i=0; i<c.c.size(); i++ ) 
    //if (c.c[i]->dim()==2)
    o << *c.c[i] << " --------------------------------------------- " << endl;
  return o;
}

/* ------------------------------------------ */

int complexs::cells_dead_or_alive ( ) const
{
  return c.size();
}

/* ------------------------------------------ */
// first 2d, then 1d, finally 0d. 
int complexs::cells (int d, int*& map1, int*& map2) const {
  int res = 0;  
  int* mp1 = new int[c.size()];
  for (int i=0; i<c.size(); i++)  mp1[i] = -1;
  for (int ii=0; ii<c.size(); ii++) {
    cell *i = c[ii];
    if (/**i &&*/ i->dim()==d) {
	      //if (*map1) 
      mp1[ii] = res;
      res++;
    } else {
      //if (*map1) (*map1)[ii] = -1;
      mp1[ii] = -1;
      }
  }
  //if (map2) {
  map2 = new int[res];
  for (int i=0; i<res; i++) map2[i] = -1;
  for (int i=0; i<c.size(); i++) 
    if (mp1[i]>=0) map2[mp1[i]] = i;
  cout<<"these are "<<d<<"th dimension cells"<<endl;
  for (int i=0; i<res; i++) {
    cout<<i<<"  "<<map2[i]<<endl;
    assert (map2[i]!=-1);
  }
  //}
  // shrink map1
  map1 = new int[res];
  for (int i=0; i<c.size(); ++i) 
    if (mp1[i]!=-1) map1[mp1[i]] = mp1[i];
  delete [] mp1;
  cout<<"number of cells="<<c.size()<<" res="<<res<<endl; 
  return res;
}

/* ------------------------------------------ */

int cell::dim() const
{
  return d;
}

/* ------------------------------------------ */

cell::operator bool() const
{
  return isalive;
}

/* ------------------------------------------ */

void complexs::remove_dead()
{
  int i=0;
  // 3 arrays: 0d, 1d and 2d
  int *newix = new int[c.size()];
  int ctr[3]={0,0,0};
  for ( int i=0; i<c.size(); i++ )
    if (*(c[i]))
      newix[i] = ctr[c[i]->d]++;
    else
      newix[i] = -1;

  std::vector<cell*> newc;
  for ( int i=0; i<c.size(); i++ )
    if (newix[i]!=-1)
      {
	newc.push_back(c[i]);
	c[i]->ID = newix[i];
	assert(newc.size()==newix[i]+1);
      }

  c = newc;

  delete[] newix;
}

/* ------------------------------------------ */

/* ------------------------------------------ */
