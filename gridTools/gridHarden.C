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
  if ( argc != 6 ) {
    printf("Set values of the potential greater than vThreshold to vNew and blur.\n");
    printf("Usage: %s srcGrid vThreshold vNew blurSteps outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  double v0 = strtod(argv[2], NULL);
  double vNew = strtod(argv[3], NULL);
  int blurSteps = atoi(argv[4]);
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded %s.\n", argv[1]);

  // Set those values greater than v0 to vNew.
  const int n = src.length();
  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    if (v > v0) src.setValue(i, vNew);
  }

  // Run some blur steps.
  for (int j = 0; j < blurSteps; j++) src.blurBigger(v0);
  
  char comments[256];
  sprintf(comments, "%s hardened to %g", inFile, vNew);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  return 0;
}
