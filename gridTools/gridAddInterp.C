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
  if ( argc < 4 ) {
    printf("Add grids of different sizes.The result has the size of the first grid.\n");
    printf("Usage: %s srcGrid0 srcGrid1... outFile\n", argv[0]);
    return 0;
  }

  const char* outFile = argv[argc-1];
  
  const int n = argc-1;
  printf("Adding %d grids.\n", n-1);

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded `%s'.\n", argv[1]);
  
  for (int i = 2; i < n; i++) {
    printf("Adding `%s'.\n", argv[i]);
    Grid g(argv[i]);
    src.addInterpolateImages(g);
  }

  char comments[256];
  sprintf(comments, "%s sum", argv[1]);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
