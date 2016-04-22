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
    printf("Usage: %s srcGrid maskGrid\n", argv[0]);
    return 0;
  }

  const char* srcGrid = argv[1];
  const char* maskGrid = argv[2];

  Grid src(srcGrid);
  Grid mask(maskGrid);

  double val = src.averageRegion(mask);
  printf("%.10g\n", val);

  return 0;
}
