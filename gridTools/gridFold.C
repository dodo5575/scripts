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
  if ( argc < 3 ) {
    printf("Fold the grid around the x axis.\n");
    printf("Usage: %s srcGrid outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const char* outFile = argv[argc-1];
  
  // Load the first grid.
  Grid src(inFile);
  printf("Loaded `%s'.\n", argv[1]);


  // Make the new grid.
  Grid out(src);
  out.zero();

  // Get the value at each node using the transformation.
  const int n = out.getSize();
  Vector3 o = out.getOrigin() - 0.5*out.getCellDiagonal();
  double len = out.getExtent().y;
  for (int i = 0; i < n; i++) {
    Vector3 r = out.getPosition(i);
    
    // Perform the transformation.
    double dy = r.y - o.y;
    double dz = r.z - o.z;
    double s = sqrt(dy*dy + dz*dz);
    double t = atan(dy/dz);
    double w = r.x - o.x;

    Vector3 r0(w, 2.0*len*t/M_PI, s);
    if (src.inGrid(r0 + o)) {
      double v0 = src.interpolatePotential(r0 + o);
      out.setValue(i, v0);
    }
  }

  char comments[256];
  sprintf(comments, "%s folded", inFile);
  out.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
