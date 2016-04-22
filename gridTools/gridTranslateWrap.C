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
  if ( argc != 6 ) {
      fprintf(stderr, "Usage: %s srcGrid dx dy dz outFile\n", argv[0]);
      fprintf(stderr, "  Translates grid boundary by (dx, dy, dz) cells, but keeps the grid\n");
      fprintf(stderr, "  values at a given coordinate the same. Assumes periodic grid.\n");
      return -1;
  }

  const char* inFile = argv[1];
  const int dx = atoi(argv[2]);
  const int dy = atoi(argv[3]);
  const int dz = atoi(argv[4]);
  const char* outFile = argv[argc-1];
  
  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);
  
  // Get grid dimensions
  int nx = src.getNx();
  int ny = src.getNy();
  int nz = src.getNz();
  
  // Make destination grid
  Grid dst(src);
  
  // Translate grid origin
  dst.setOrigin(src.getOrigin() + src.getBasis().transform(Vector3(dx, dy, dz)));
  
  // Translate grid values
  for (int sx = 0; sx < nx; sx++) {
      for (int sy = 0; sy < ny; sy++) {
	  for (int sz = 0; sz < nz; sz++) {
	      // calculate destination indices (x y z) from source indices (sx sy sz)
	      // subtract because (dx dy dz) indicates movement of the grid _boundary_,
	      // so relative to the boundary, the grid values move the opposite way
	      int x = (sx + nx - dx) % nx;
	      int y = (sy + ny - dy) % ny;
	      int z = (sz + nz - dz) % nz;
	      double val = src.getValue(sx, sy, sz);
	      dst.setValue(x, y, z, val);
	  }
      }
  }
  
  char comments[256];
  sprintf(comments, "(%s) translated by %d %d %d", inFile, dx, dy, dz);
  dst.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
