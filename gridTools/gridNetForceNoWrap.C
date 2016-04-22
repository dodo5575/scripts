// Author: Jeff Comer <jcomer2@illinois.edu>
// Apply a uniform external field to a grid.
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
  if ( argc != 4 ) {
    printf("Compute the net force for a set of point charges in the grid.\n");
    printf("Usage: %s srcGrid chargeFile coordinateFile\n", argv[0]);
    return 0;
  }

  Grid src(argv[1]);
  Scatter charge(1, argv[2]);
  Scatter pos(argv[3]);

  const int n = pos.length();
  const int nc = charge.length();
  if (n != nc) {
    printf("Error! chargeFile and coordinateFile have different lengths.\n");
    return 0;
  }

  Vector3 sumForce = Vector3(0.0);
  double sumPot = 0.0;
  double v;
  Vector3 f;
  for (int i = 0; i < n; i++) {
    double q = charge.get(i).x;
    Vector3 r = pos.get(i);
    if (src.inGridInterp(r)) {
      src.interpolate(&v, &f, r);
      sumPot += q*v;
      sumForce += q*f;
    }
  }

  printf("%.12g %.12g %.12g %.12g\n", sumPot, sumForce.x, sumForce.y, sumForce.z);

  return 0;
}
