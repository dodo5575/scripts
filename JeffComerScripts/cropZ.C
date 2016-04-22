#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
//#include "Scatter.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Usage: %s srcFile z0 z1 outputFile\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  double z0 = strtod(argv[2],NULL);
  double z1 = strtod(argv[3],NULL);

  // Find the z indices.
  Vector3 l0 = src.transformTo(Vector3(0.0,0.0,z0));
  Vector3 l1 = src.transformTo(Vector3(0.0,0.0,z1));
  int iz0 = int(floor(l0.z + 0.5));
  int iz1 = int(floor(l1.z + 0.5));
  if (iz0 > iz1) {
    int i = iz0;
    iz0 = iz1;
    iz1 = i;
  }

  // Crop it.
  src.crop(0, 0, iz0, src.getNx()-1, src.getNy()-1, iz1); 

  src.write(argv[argc-1], "touched up");
  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
