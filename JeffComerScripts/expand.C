#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s srcFile factor outputFile\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  double factor = strtod(argv[2],NULL);

  int nx0 = src.getNx();
  int ny0 = src.getNy();
  int nz0 = src.getNz();

  int nx = int(factor*nx0);
  int ny = int(factor*ny0);
  int nz = int(factor*nz0);
  Vector3 origin = factor*src.getOrigin();
  Grid dst(src.getBasis(), origin, nx, ny, nz);

  int sx0 = (nx - nx0)/2;
  int sy0 = (ny - ny0)/2;
  int sz0 = (nz - nz0)/2;

  for (int ix = 0; ix < nx; ix++) {
      for (int iy = 0; iy < ny; iy++) {
	  for (int iz = 0; iz < nz; iz++) {
	    int jx = ix - sx0;
	    int jy = iy - sy0;
	    int jz = iz - sz0;

	    if (jx < 0) jx = 0;
	    if (jx >= nx0) jx = nx0-1;
	    if (jy < 0) jy = 0;
	    if (jy >= ny0) jy = ny0-1;
	    if (jz < 0) jz = 0;
	    if (jz >= nz0) jz = nz0-1;

	    double v = src.getValue(jx,jy,jz);
	    dst.setValue(ix,iy,iz,v);
	  }
      }
  }

  char comments[256];
  sprintf(comments, "%s expanded by %g", argv[1], factor);

  dst.write(argv[argc-1], comments);
  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
