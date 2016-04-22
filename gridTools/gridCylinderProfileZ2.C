// Author: David Wells <dbwells2@illinois.edu>
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
  if ( argc != 6 ) {
    printf("Extract a radial profile of the grid, averaged in angle and z and centered on (xo,yo).\n");
    printf("Usage: %s srcGrid xo yo radius outFile\n", argv[0]);
    return 0;
  }
  
  const char* inFile = argv[1];
  const double xo = strtod(argv[2], NULL);
  const double yo = strtod(argv[3], NULL);
  const double radius = strtod(argv[4], NULL);
  const char* outFile = argv[argc-1];
  
  // Load the first grid.
  Grid src(inFile);
  printf("Loaded `%s'.\n", inFile); 
  
  src.cylinderProfileZ2(outFile, xo, yo, radius);
  return 0;
}
