#ifndef _SEGMENTUNION_H
#define _SEGMENTUNION_H

#include <vector>
using namespace std;

namespace util {
class Seg {
  public:
  double s, e;
  Seg() {}
  Seg(double s1, double e1) :s(s1), e(e1) {} 
};

class SegmentUnion {
  public:
  vector<Seg> segV;

  SegmentUnion() {}
  void addPair(double s, double e) {
    Seg seg(s,e);
    segV.push_back(seg);
  }
  void addPair(Seg& seg) {
    addPair(seg.s,seg.e);
  }
  int nbSeg() {
    return segV.size();
  }

  void adjust(int l, int r) {
    int i, j;
    double x, y;
    if (l < r){
       x = segV[l].s; y = segV[l].e;
       i = l;
       j = r;    
       while(i<j) {
         while (i<j && segV[j].s>=x) j--;
         if(i<j) { 
           segV[i].s   = segV[j].s;   
           segV[i++].e = segV[j].e;
         }
         while (i<j && segV[i].s<=x) i++;
         if(i<j) {
           segV[j].s   = segV[i].s;       
           segV[j--].e = segV[i].e;
         }
       }
       segV[i].s = x; segV[i].e = y;        
       adjust(l,i-1);
       adjust(i+1,r);
    }
  }
  
  void unite() {
    int n = segV.size();
    if (n==0) return;
    adjust(0, n-1);
    vector<Seg> temp;
    double x1 = segV[0].s, x2 = segV[0].e;
    if (n==1) temp.push_back(segV[0]);
    else {
      for (int i=1; i<n; ++i) {
        if(segV[i].s>x2) { 
          Seg seg(x1,x2);
          temp.push_back(seg);
          x1 = segV[i].s; x2 = segV[i].e;
        }   
        else if(segV[i].e>x2) x2 = segV[i].e;
      }
    }
    Seg segf(x1,x2);
    temp.push_back(segf);
    segV.swap(temp);
  }
};
}


#endif
