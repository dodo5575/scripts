// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 5 && argc != 3) {
    printf("Get the force at a particular point or points in the potential grid.\n");
    printf("Usage: %s srcGrid x y z\n", argv[0]);
    printf("Usage: %s srcGrid inputPositionFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  // Load the grid.
  Grid src(inFile);

  Scatter* posList;

   if (argc == 5) {
    // Just do one point.
    double x = strtod(argv[2], NULL);
    double y = strtod(argv[3], NULL);
    double z = strtod(argv[4], NULL);
    Vector3 r0(x, y, z);

    posList = new Scatter(&r0, 1);
  } else {
    // Load the points from the file.
    posList = new Scatter(argv[2]);
  }

  const int n = posList->length();
  double v;
  Vector3 f;
  for (int i = 0; i < n; i++) { 
    src.interpolate(&v, &f, posList->get(i));
    printf("%.12g %.12g %.12g %.12g\n", v, f.x, f.y, f.z);
  }

  delete posList;

  return 0;
}
