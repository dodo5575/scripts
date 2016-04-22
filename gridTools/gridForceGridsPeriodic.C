// Author: Jeff Comer <jcomer2@illinois.edu>
// Generate force grids with periodic images on all sides.
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
    printf("Generate force grids with periodic images on all sides.\n");
    printf("Usage: %s srcGrid outPrefix\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const char* outPrefix = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);
  Grid big(src.tile(3, 3, 3));

  // Compute the forces.
  Grid* fx = new Grid(big);
  Grid* fy = new Grid(big);
  Grid* fz = new Grid(big);
  Grid* fm = new Grid(big);

  //Grid* v = new Grid(src);
  //Grid* fx = NULL;
  //Grid* fy = NULL;
  //Grid* fz = NULL;
  //Grid* fm = NULL;

    //printf("pointer: %d\n", fx);
  printf("Computing force...\n");
  src.computeForce(fx, fy, fz, fm);
  //printf("pointer: %d\n", fx);

  char comments[256];
  sprintf(comments, "%s force", inFile);
  printf("%s.\n", comments);
  printf("Writing...\n");

  // Write the force grids.
  char outX[256];
  sprintf(outX, "%s.fx.dx", outPrefix);
  //fx->write(outX, comments);

  char outY[256];
  sprintf(outY, "%s.fy.dx", outPrefix);
  //fy->write(outY, comments);

  char outZ[256];
  sprintf(outZ, "%s.fz.dx", outPrefix);
  fz->write(outZ, comments);
  
  char outM[256];
  sprintf(outM, "%s.fm.dx", outPrefix);
  //fm->write(outM, comments);

  // Make profiles.
  sprintf(outZ, "%s.fz.dat", outPrefix);
  fz->profileZ(0.5, 0.5, outZ);

  sprintf(outM, "%s.fm.dat", outPrefix);
  //fm->profileZ(0.5, 0.5, outM);

  //char outV[256];
  //sprintf(outV, "%s.v.dx", outPrefix);
  //v->write(outV, comments);

  printf("Wrote output files.\n");

  delete fx;
  delete fy; 
  delete fz; 
  delete fm; 
  //delete v;
  return 0;
}
