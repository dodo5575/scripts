// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

double min(double x, double y) {
  if (x <= y) return x;
  return y;
}

double max(double x, double y) {
  if (x >= y) return x;
  return y;
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Shift the potential so that the average of the region z0 <= z < z1 is zero.\n");
    printf("Usage: %s srcGrid0 z0 z1 outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const double z0 = strtod(argv[2], NULL);
  const double z1 = strtod(argv[3], NULL);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded `%s'.\n", inFile); 

  // Get the region.
  Vector3 diag = src.getCellDiagonal();
  Vector3 p0 = src.getOrigin() - diag;
  Vector3 p1 = src.getDestination() + diag;

  Vector3 r0;
  r0.x = min(p0.x, p1.x);
  r0.y = min(p0.y, p1.y);
  r0.z = z0;
  
  Vector3 r1;
  r1.x = max(p0.x, p1.x);
  r1.y = max(p0.y, p1.y);
  r1.z = z1;

  printf("Region: %g <= z < %g\n", r0.z, r1.z);

  // Get the average.
  double v0 = src.averageRegion(r0, r1);
  printf("Average value: %g\n", v0);
  
  // Shift the potential.
  src.shift(-v0);   

  char comments[256];
  sprintf(comments, "%s shifted", inFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
