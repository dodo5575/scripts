// Author: Jeff Comer <jcomer2@illinois.edu>
// Sets the membrane to 1, everything else to zero.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

// Is a point in the material making up the membrane? 
bool inPoreWalls(Vector3 r, double l, double sx, double sy, double slope, double paraZ) {
  if (fabs(r.z) > l) return false;

  if (fabs(r.z) < paraZ) {
    // parabolic region
     double dx = r.x/(sx + 0.5*slope/paraZ*r.z*r.z);
     double dy = r.y/(sy + 0.5*slope/paraZ*r.z*r.z);
     if (dx*dx + dy*dy < 1) return false;
  } else {
    // linear region
    double dx = r.x/(sx - 0.5*slope*paraZ + slope*fabs(r.z));
    double dy = r.y/(sy - 0.5*slope*paraZ + slope*fabs(r.z));
    if (dx*dx + dy*dy < 1) return false;
  }
  
  return true;
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc < 8 ) {
    printf("Usage: %s inGrid poreLength poreDiameterX poreDiameterY poreAngle parabolicLen outGrid\n", argv[0]);
    printf("Make a grid with the dimensions of inGrid containing a double cone pore.\n");
    return 0;
  }
  
  const char* inGrid = argv[1];
  double poreLength = strtod(argv[2], NULL);
  double poreDiameterX = strtod(argv[3], NULL);
  double poreDiameterY = strtod(argv[4], NULL);
  double poreAngle = strtod(argv[5], NULL);
  double parabolicLen = strtod(argv[6], NULL);
  const char* outGrid = argv[argc-1];

  // Make the grid.
  Grid src(inGrid);

  // Get the geometry.
  double l = 0.5*poreLength;
  double sx = 0.5*poreDiameterX;
  double sy = 0.5*poreDiameterY;
  double pi = 4.0*atan(1.0);
  double slope = tan(pi*poreAngle/180.0);
  double paraZ = 0.5*parabolicLen;

  printf("Scanning points\n");
  int count = 0;
  const int n = src.length();
  for (int i = 0; i < n; i++) {
    Vector3 r(src.getPosition(i));
    
    if (inPoreWalls(r, l, sx, sy, slope, paraZ)) {
      src.setValue(i, 1.0);
      count++;
    }
  }
 
  double totalVol = src.getVolume();
  double srcVol = (totalVol*count)/n;
  double remainVol = totalVol - srcVol;
  printf("Total points: %d\n", n);
  printf("Membrane points: %d\n", count);
  printf("Remaining points: %d\n", n-count);
  printf("Membrane fraction: %.10g\n", double(count)/n);
  printf("Total volume: %.10g\n", totalVol);
  printf("Membrane volume: %.10g\n", srcVol);
  printf("Remaining volume: %.10g\n", remainVol);
  
  char comments[256];
  sprintf(comments, "poreLength %g poreDiameterX %g poreDiameterY %g poreAngle %g parabolicLen %g",
	  poreLength, poreDiameterX, poreDiameterY, poreAngle, parabolicLen);
  src.write(outGrid, comments);
  printf("%s\n", comments);
  printf("Wrote `%s'.\n", outGrid);
  
  
  return 0;
}
