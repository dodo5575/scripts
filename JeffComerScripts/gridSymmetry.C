// Author: Jeff Comer <jcomer2@illinois.edu>
// Average a grid around the z-axis a given number of times.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc < 4 ) {
    printf("Usage: %s inGrid symmetry outGrid\n", argv[0]);
    return 0;
  }
  
  Grid grid(argv[1]);
  Grid gridOut(grid);
  Grid gridIn(grid);
  int symmetry = atoi(argv[2]);

  const double twoPi = 8.0*atan(1.0);
  const int n = grid.length();
  Matrix3 box = grid.getBox();
  double lx = box.ex().length();
  double ly = box.ey().length();
  double rad = (lx < ly)?lx:ly;

  // Make gridOut the potential for the region outside the cylinder.
  // Make gridIn the potential for the region inside the cylinder.
  for (int i = 0; i < n; i++) {
    Vector3 r = grid.getPosition(i);
    if (r.length() < rad) gridOut.setValue(i, 0.0);
    else gridIn.setValue(i, 0.0);
  }

  // Make the rotation matrix.
  Grid gridSum(gridIn);
  Vector3 ex(cos(twoPi/symmetry), sin(twoPi/symmetry), 0.0);
  Vector3 ey(-sin(twoPi/symmetry), cos(twoPi/symmetry), 0.0);
  Vector3 ez(0.0, 0.0, 1.0);
  Matrix3 rotate(ex, ey, ez);

  // Rotate gridIn around the symmetry axes and add it to the sum.
  for (int j = 1; j < symmetry; j++) {   
    gridIn.setBasis(rotate*gridIn.getBasis());
    gridIn.setOrigin(rotate.transform(gridIn.getOrigin()));
    gridSum.addInterpolate(gridIn);
  }
  gridSum.scale(1.0/symmetry);
  gridSum.add(gridOut);

  char comments[256];
  sprintf(comments, "%s with %d symmetry", argv[1], symmetry);
  gridSum.write(argv[argc-1], comments);

  return 0;
}
