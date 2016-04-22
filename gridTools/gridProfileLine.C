// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 10 ) {
    printf("Extract a profile of the grid along (x0,y0,z0) and (x1,y1,z1) consisting of n+1 points.\n");
    printf("Usage: %s inGrid n x0 y0 z0 x1 y1 z1 outFile\n", argv[0]);
    return 0;
  }

  const char* inGrid = argv[1];
  const int n = atoi(argv[2]);
  const double x0= strtod(argv[3], NULL);
  const double y0 = strtod(argv[4], NULL);
  const double z0 = strtod(argv[5], NULL);
  const double x1 = strtod(argv[6], NULL);
  const double y1 = strtod(argv[7], NULL);
  const double z1 = strtod(argv[8], NULL);
  const char* outFile = argv[argc-1];


  // Get the geometry.
  Vector3 r0(x0,y0,z0);
  Vector3 r1(x1,y1,z1);
  Vector3 dr = (r1-r0)/n;
  double dist = dr.length();

  // Load the first grid.
  Grid src(inGrid);
  printf("Loaded `%s'.\n", inGrid);
  printf("Generating profile of %d+1 points along:\n", n);
  printf("%s\n%s\n", r0.toString().val(), r1.toString().val());

  // Write the profile.
  FILE* out = fopen(outFile, "w");
  for (int i = 0; i <= n; i++) {
    Vector3 r = r0 + dr*i;
    double v = src.interpolatePotential(r);

    fprintf(out, "%.10g %.10g\n", i*dist, v);
  }
  fclose(out);
  
  return 0;
}
