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
    printf("  Prints the grid values for the cell at the given (integer) coordinate\n");
    printf("  Assumes the grid is periodic in all dimensions\n");
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
  
  for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
	  for (int k = 0; k < 2; k++) {
	      // get grid value at coordinate
	      int wx = (x + i) % nx;
	      int wy = (y + j) % ny;
	      int wz = (z + k) % nz;
	      double val = src.getValue(wx, wy, wz);
	      printf("(%d %d %d) : %f\n", wx, wy, wz, val);
	  }
      }
  }
  return 0;
}
