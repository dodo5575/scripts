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
  if ( argc != 3 ) {
    printf("out = 1 / out\n");
    printf("Usage: %s srcGrid outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const char* outFile = argv[argc-1];

  // Load the grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  src.invert();
  
  char comments[256];
  sprintf(comments, "1 / (%s)", inFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
