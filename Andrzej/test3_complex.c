
#include <iostream>
#include <complex.h>

using namespace std;


/* ------------------------------------------ */

int main ()
{
  complex c("t.t");

  cout << c << endl;// << c.check() << endl;

  cout << " ------------------------------------------------------- " << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << " ------------------------------------------------------- " << endl;

  c.get_cell(0)->merge(c.get_cell(1));
  //c.remove_dead();
  cout << c << endl;
  /*cout << c.check() << endl;

  cout << " ------------------------------------------------------- " << endl;

  c.get_cell(0)->merge(c.get_cell(3));
  cout << c << endl;
  cout << c.check() << endl;

  c.remove_dead();

  cout << " ------------------------------------------------------- " << endl;
  cout << c << endl;
  cout << c.check() << endl;*/

  return 0;
}
