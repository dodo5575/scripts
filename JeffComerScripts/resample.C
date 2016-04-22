#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainResample(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s inGrid factor outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  double factor = strtod(argv[2],NULL);

  int nx = int(floor(orig.getNx()*factor));
  int ny = int(floor(orig.getNy()*factor));
  int nz = int(floor(orig.getNz()*factor));

  Grid resamp = orig.resample(nx, ny, nz);

  char comments[256];
  sprintf(comments, "%s resampled by %.10g", argv[1], factor);
  printf("%s\n", comments);
  resamp.write(argv[argc-1], comments);
  printf("Wrote %s.\n", argv[argc-1]);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainResample(argc, argv);
}
