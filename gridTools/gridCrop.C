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
  if ( argc != 9 ) {
    printf("Crop the grid.\n");
    printf("Usage: %s srcGrid x0 y0 z0 x1 y1 z1 outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  double x0 = strtod(argv[2],NULL);
  double y0 = strtod(argv[3],NULL);
  double z0 = strtod(argv[4],NULL);
  double x1 = strtod(argv[5],NULL);
  double y1 = strtod(argv[6],NULL);
  double z1 = strtod(argv[7],NULL);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);
  
  src.crop(Vector3(x0,y0,z0), Vector3(x1,y1,z1));

  char comments[256];
  sprintf(comments, "%s cropped", inFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
