// Author: Jeff Comer <jcomer2@illinois.edu>
// Apply a uniform external field to a grid.
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Usage: %s srcGrid nx ny nz outGrid\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const int nx = atoi(argv[2]);
  const int ny = atoi(argv[3]);
  const int nz = atoi(argv[4]);
  const char* outName = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  // Compute the forces.
  Grid inter(src, nx, ny, nz);
  
  char comments[256];
  sprintf(comments, "%s resamp", inFile);
  printf("%s.\n", comments);
  
  // Write the result.
  inter.write(outName, comments);
  printf("Wrote the output file `%s'.\n", outName);

  return 0;
}
