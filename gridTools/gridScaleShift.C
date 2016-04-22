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
  if ( argc != 5 ) {
    printf("out = scale*src + shift\n");
    printf("Usage: %s srcGrid scale shift outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const double scale = strtod(argv[2],NULL);
  const double shift = strtod(argv[3],NULL);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  src.scale(scale);
  src.shift(shift);
  
  char comments[256];
  sprintf(comments, "(%s)*%g + %g", inFile, scale, shift);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
