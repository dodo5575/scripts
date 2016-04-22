// Author: Jeff Comer <jcomer2@illinois.edu>
// Author: David Wells <dbwells2@illinois.edu>
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
    printf("Usage: %s srcGrid x y z\n", argv[0]);
    printf("  Prints the grid value at the given (integer) coordinate\n");
    return 0;
  }

  const char* srcGrid = argv[1];
  const int x = atoi(argv[2]);
  const int y = atoi(argv[3]);
  const int z = atoi(argv[4]);
  Grid src(srcGrid);
  
  // get grid dimensions
  int nx = src.getNx();
  int ny = src.getNy();
  int nz = src.getNz();
  
  fprintf(stderr, "grid dimensions: %d %d %d\n", nx, ny, nz);
  
  // get grid value at coordinate
  double val = src.getValue(x, y, z);
  
  // print
  if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz) {
      printf("%f\n", val);
      return 0;
  } else {
      fprintf(stderr, "Coordinate not within grid!\n");
      return -1;
  }
}
