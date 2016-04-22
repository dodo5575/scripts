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
  if ( argc != 3 ) {
    printf("Write the grid nodes as carbon atoms in the VMD readable XYZ format.\n");
    printf("Usage: %s srcGrid outFile\n", argv[0]);
    return 0;
  }

  const char* srcGrid = argv[1];
  const char* outFile = argv[argc-1];
  
  Grid src(srcGrid);
  const int n = src.length();

  // Open the file and write the header.
  FILE* out = fopen(outFile, "w");
  fprintf(out, "%d\n", src.length());
  fprintf(out, " xyz file\n");

  // Write the nodes.
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);
    fprintf(out, " C %.10g %.10g %.10g\n", r.x, r.y, r.z);
  }
  fclose(out);
  

  return 0;
}
