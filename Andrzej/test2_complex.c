
#include <iostream>
#include <complex.h>

using namespace std;

/* ------------------------------------------ */

int main ( )
{
  complex c;

  int c0[8];
  int c1[12];
  int c2[6];

  for ( int i=0; i<8; i++ )
    c0[i] = c.add_cell ( new cell(0,0) );
  for ( int i=0; i<12; i++ )
    c1[i] = c.add_cell ( new cell(1,1) );
  for ( int i=0; i<6; i++ )
    c2[i] = c.add_cell ( new cell(2,1) );

  c.makeface(c1[0],c0[0],-1);
  c.makeface(c1[0],c0[1],1);
  c.makeface(c1[1],c0[1],-1);
  c.makeface(c1[1],c0[2],1);
  c.makeface(c1[2],c0[2],-1);
  c.makeface(c1[2],c0[3],1);
  c.makeface(c1[3],c0[0],-1);
  c.makeface(c1[3],c0[3],1);
  c.makeface(c1[4],c0[4],-1);
  c.makeface(c1[4],c0[5],1);
  c.makeface(c1[5],c0[5],-1);
  c.makeface(c1[5],c0[6],1);
  c.makeface(c1[6],c0[6],-1);
  c.makeface(c1[6],c0[7],1);
  c.makeface(c1[7],c0[7],-1);
  c.makeface(c1[7],c0[4],1);
  c.makeface(c1[8],c0[0],-1);
  c.makeface(c1[8],c0[4],1);
  c.makeface(c1[9],c0[1],-1);
  c.makeface(c1[9],c0[5],1);
  c.makeface(c1[10],c0[2],-1);
  c.makeface(c1[10],c0[6],1);
  c.makeface(c1[11],c0[3],-1);
  c.makeface(c1[11],c0[7],1);

  c.makeface(c2[0],c1[3],1);
  c.makeface(c2[0],c1[8],-1);
  c.makeface(c2[0],c1[7],1);
  c.makeface(c2[0],c1[11],1);
  c.makeface(c2[1],c1[0],1);
  c.makeface(c2[1],c1[1],1);
  c.makeface(c2[1],c1[2],1);
  c.makeface(c2[1],c1[3],-1);
  c.makeface(c2[2],c1[2],1);
  c.makeface(c2[2],c1[6],-1);
  c.makeface(c2[2],c1[10],-1);
  c.makeface(c2[2],c1[11],1);
  c.makeface(c2[3],c1[0],1);
  c.makeface(c2[3],c1[4],-1);
  c.makeface(c2[3],c1[8],-1);
  c.makeface(c2[3],c1[9],1);
  c.makeface(c2[4],c1[9],-1);
  c.makeface(c2[4],c1[1],1);
  c.makeface(c2[4],c1[5],-1);
  c.makeface(c2[4],c1[10],1);
  c.makeface(c2[5],c1[4],1);
  c.makeface(c2[5],c1[5],1);
  c.makeface(c2[5],c1[6],1);
  c.makeface(c2[5],c1[7],1);

  c.get_cell(c2[0])->merge(c.get_cell(c2[2]));
  cout << c << endl;
  cout << c.check() << endl;

  cout << " ------------------------------------------------------- " << endl;

  c.get_cell(c2[0])->merge(c.get_cell(c2[2]));
  c.get_cell(c2[0])->merge(c.get_cell(c2[3]));
  c.get_cell(c2[0])->merge(c.get_cell(c2[4]));

  cout << " ------------------------------------------------------- " << endl;

  cout << c << endl;
  cout << c.check() << endl;

  c.remove_dead();

  cout << " ------------------------------------------------------- " << endl;
  cout << c << endl;
  cout << c.check() << endl;

  return 0;
}
