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
  if ( argc != 3  ) {
    printf("Subtract grids of different sizes (srcGrid0 - srcGrid1). The result has the size of the first grid.\n");
    printf("Usage: %s srcGrid0 srcGrid1 outFile\n", argv[0]);
    return 0;
  }

  const char* outFile = argv[argc-1];
  
  const int n = argc-1;
  printf("Adding %d grids.\n", n-1);

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded `%s'.\n", argv[1]);
  
  // Load the second grid.
  Grid g(argv[2]);
  printf("Loaded `%s'.\n", argv[2]);
  g.scale(-1.0);
  src.addInterpolateImages(g);

  char comments[256];
  snprintf(comments, 256, "%s sum", argv[1]);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
