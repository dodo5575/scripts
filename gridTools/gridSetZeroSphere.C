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
  if ( argc != 4 ) {
    printf("Shift the potential so that the average of the region r < rad is zero.\n");
    printf("Usage: %s srcGrid0 rad outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const double rad = strtod(argv[2], NULL);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded `%s'.\n", inFile); 

  printf("Region: r < %g\n", rad);

  // Get the average.
  double v0 = src.averageSphere(Vector3(0.0), rad);
  printf("Average value: %g\n", v0);
  
  // Shift the potential.
  src.shift(-v0);   

  char comments[256];
  sprintf(comments, "%s shifted", inFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
