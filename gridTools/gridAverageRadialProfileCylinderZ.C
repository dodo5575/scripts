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
  if ( argc != 5 ) {
    printf("Extract a radial profile of the grid, averaged in angle and z.\n");
    printf("Usage: %s srcGrid nbins rmax outFile\n", argv[0]);
    return 0;
  }
  
  const char* inFile = argv[1];
  int nbins = strtol(argv[2], NULL, 0);
  double rmax = strtod(argv[3], NULL);
  const char* outFile = argv[argc-1];
  
  // Load the first grid.
  Grid src(inFile);
  printf("Loaded `%s'.\n", inFile); 
  
  src.averageRadialProfileCylinderZ(outFile, nbins, rmax);
  
  return 0;
}
