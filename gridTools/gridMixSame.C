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
  if ( argc < 5 || argc % 2 != 1) {
    printf("Compose grids of the same size using masks.\n");
    printf("Usage: %s baseGrid srcGrid0 maskGrid0 srcGrid1 maskGrid1... outFile\n", argv[0]);
    return 0;
  }

  const char* outFile = argv[argc-1];
  
  const int n = argc-1;
  printf("Composing %d grids.\n", (n-1)/2);

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded `%s'.\n", argv[1]);
  
  for (int i = 2; i < n; i+=2) {
    printf("Mixing `%s' with mask `%s'.\n", argv[i], argv[i+1]);
    Grid g(argv[i]);
    Grid m(argv[i+1]);
    bool good = src.mix(g, m);

    if (!good) printf("Warning! `%s' is not the same size as `%s'. Grid not added.\n", argv[i], argv[1]);  }

  char comments[256];
  snprintf(comments, 256, "%s sum", argv[1]);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
