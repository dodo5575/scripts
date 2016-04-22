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
  if ( argc != 4 ) {
    printf("Usage: %s inGrid0 inGrid1 outGrid\n", argv[0]);
    return 0;
  }

  const char* inGrid0 = argv[1];
  const char* inGrid1 = argv[2];
  const char* outGrid = argv[argc-1];

  // Load the first grid.
  Grid src0(inGrid0);
  printf("Loaded %s.\n", inGrid0);
  Grid src1(inGrid1);
  printf("Loaded %s.\n", inGrid1);

  Grid dest(src0);

  const int n = src0.length();
  for (int i = 0; i < n; i++) {
    Vector3 r = src0.getPosition(i);
    
    Vector3 f0 = src0.interpolateForce(r);
    Vector3 f1 = src1.interpolateForce(r);
    
    dest.setValue(i, f0.dot(f1));
  }
  
  char comments[256];
  snprintf(comments, 256, "force(%s) dot force(%s)", inGrid0, inGrid1);
  dest.write(outGrid, comments);
  printf("Wrote `%s'.\n", outGrid);

  return 0;
}
