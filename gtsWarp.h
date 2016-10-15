#ifndef _GTSWARP_H
#define _GTSWARP_H

#include "gts.h"

bool isSame(GtsEdge* e1, GtsEdge* e2) {
  GtsVertex *v11 = GTS_SEGMENT(e1)->v1;
  GtsVertex *v12 = GTS_SEGMENT(e1)->v2;
  GtsVertex *v21 = GTS_SEGMENT(e2)->v1;
  GtsVertex *v22 = GTS_SEGMENT(e2)->v2;
  if (v11==v21 && v12==v22 || v11==v22 && v12==v21) return true;
  return false;
}

GtsSurface* toGts(Mesh& m) {
  GtsSurface *surface = gts_surface_new(gts_surface_class(),gts_face_class(),
    gts_edge_class(),gts_vertex_class());
  // first add vertices
  int i, j, nv = m.NbVertex();
  GtsVertex **v = new GtsVertex*[nv];
  for (i=0; i<nv; ++i) 
    v[i] = gts_vertex_new(gts_vertex_class(),
      m.v(i).x(),m.v(i).y(),m.v(i).z());
  // then add edges and faces
  int nf = m.NbFace();
  GtsFace **f = new GtsFace*[nf];
  vector<GtsEdge*> ea;
  // 3 temp edges for creating a face
  GtsEdge **ef = new GtsEdge*[3];
  for (i=0; i<nf; ++i) {
    // traverse edges
    for (j=0; j<3; ++j) {
      GtsEdge* e = gts_edge_new(gts_edge_class(), 
        v[m.f(i).v(j)], v[m.f(i).v((j+1)%3)]);
      if (ea.empty()) {ef[0] = e; ea.push_back(e); continue;}
      // check if this edge has been added or not
      int nEa = ea.size();
      bool isNew = true;
      for (int k=0; k<nEa; ++k) 
        if (isSame(ea[k],e)) {
          // destroy edge e
          gts_object_destroy(GTS_OBJECT(e));
          ef[j] = ea[k];
          isNew = false;
          break;
        }
      // add edge e into the edge array
      if (isNew) {
        ea.push_back(e);
        ef[j] = e;
      }
    }
    // create face
    f[i] = gts_face_new(gts_face_class(),ef[0],ef[1],ef[2]);
    gts_surface_add_face(surface,f[i]);
  }
  return surface;
}

static void writeVertex(GtsPoint *p, gpointer *data) {
  g_hash_table_insert((GHashTable*)data[1],p,
    GUINT_TO_POINTER(++(*((guint*)data[0]))));
  ((Mesh*)data[2])->addVertex(p->x,p->y,p->z);
}

static void writeFace(GtsTriangle *t, gpointer *data) {
  // get indices of 3 vertices in this trianlge
  GtsVertex *v1  = GTS_SEGMENT(t->e1)->v1;
  GtsVertex *v2  = GTS_SEGMENT(t->e1)->v2;
  GtsVertex *v21 = GTS_SEGMENT(t->e2)->v1;
  GtsVertex *v22 = GTS_SEGMENT(t->e2)->v2;
  GtsVertex *v3;
  if (v21!=v1 && v21!=v2) v3 = v21;
  else v3 = v22;
  int vi1 = GPOINTER_TO_UINT(g_hash_table_lookup((GHashTable*)data[1],v1))-1;
  int vi2 = GPOINTER_TO_UINT(g_hash_table_lookup((GHashTable*)data[1],v2))-1;
  int vi3 = GPOINTER_TO_UINT(g_hash_table_lookup((GHashTable*)data[1],v3))-1;
  ((Mesh*)data[2])->addFace(vi1,vi2,vi3);
}

void toMesh(GtsSurface *s, Mesh& m) {
  guint n = 0;
  GHashTable *vind;
  gpointer data[3];
  data[0] = &n;
  data[1] = vind = g_hash_table_new(NULL,NULL);
  data[2] = &m;
  // traverse vertices to build vertex table
  gts_surface_foreach_vertex(s,(GtsFunc)writeVertex,data);
  // traverse triangles to build face table
  gts_surface_foreach_face(s,(GtsFunc)writeFace,data);
  m.buildAdjacency();
  m.orientFaces();
  m.buildAdjacency();
  m.computeNormal();
}

void coarsen(Mesh& in, Mesh& out, unsigned int stopNo) {
  GtsSurface *s = toGts(in);
  guint stop_number = stopNo;
  gdouble fold = 3.1415926/180.0;
  GtsVolumeOptimizedParams params = {0.5, 0.5, 0};
  GtsStopFunc stopFunc = (GtsStopFunc) gts_coarsen_stop_number;
  GtsKeyFunc costFunc = (GtsKeyFunc)gts_volume_optimized_cost;
  GtsCoarsenFunc coarsenFunc = NULL;
    //(GtsCoarsenFunc)gts_volume_optimized_cost;
  gpointer costData = &params;
  gpointer coarsenData = NULL;//&params;
  gpointer stopData = &stop_number;
  if (gts_surface_is_manifold(s)) cout<<"Manifold"<<endl;
  else cout<<"Not Manifold"<<endl;
  gts_surface_coarsen(s,
                      costFunc, costData, 
                      coarsenFunc, coarsenData,
                      stopFunc, stopData, 
                      fold);
  if (gts_surface_is_manifold(s)) cout<<"Manifold"<<endl;
  else cout<<"Not Manifold"<<endl;
  toMesh(s,out);
  cout<<out.NbVertex()<<"  "<<out.NbFace()<<endl;
  gts_object_destroy(GTS_OBJECT(s));
}

#endif
