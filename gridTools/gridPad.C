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
  if ( argc != 9 ) {
    printf("Usage: %s srcGrid padX0 padY0 padZ0 padX1 padY1 padZ1 outGrid\n", argv[0]);
    return 0;
  }

  const char* srcGrid = argv[1];
  const int padX0 = atoi(argv[2]);
  const int padY0 = atoi(argv[3]);
  const int padZ0 = atoi(argv[4]);
  const int padX1 = atoi(argv[5]);
  const int padY1 = atoi(argv[6]);
  const int padZ1 = atoi(argv[7]);
  const char* outFile = argv[argc-1];

  Grid src(srcGrid);
  Grid dest = src.pad(padX0, padY0, padZ0, padX1, padY1, padZ1);

  char comments[256];
  sprintf(comments, "%s padded by %d %d %d", srcGrid, padX1, padY1, padZ1);
  dest.write(outFile, comments);

  return 0;
}
