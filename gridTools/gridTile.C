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
  if ( argc != 6 ) {
    printf("Usage: %s srcGrid nx ny nz outGrid\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const int nx = atoi(argv[2]);
  const int ny = atoi(argv[3]);
  const int nz = atoi(argv[4]);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", inFile);
  Grid tiled = src.tile(nx,ny,nz);
  
  char comments[256];
  sprintf(comments, "tile %s by %d %d %d", inFile, nx, ny, nz);
  tiled.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
