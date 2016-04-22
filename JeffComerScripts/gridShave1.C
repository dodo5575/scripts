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
  if ( argc != 4 ) {
    printf("Usage: %s srcGrid count outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const int count = atoi(argv[2]);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  int nx = src.getNx();
  int ny = src.getNy();
  int nz = src.getNz();
  Matrix3 basis = src.getBasis();
  Vector3 origin = src.getOrigin();
  nx-=count;
  ny-=count;

  Grid out(basis, origin, nx, ny, nz);
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int iz = 0; iz < nz; iz++) {
	double v = src.getValue(ix, iy, iz);
	out.setValue(ix, iy, iz, v);
      }
    }
  }
  
  char comments[256];
  sprintf(comments, "%s shaved", inFile);
  out.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
