#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 4 && argc != 6 ) {
    printf("Make a new grid (all zeros) from a basis file and desired grid spacing.\n");
    printf("basisFile lists the cell basis vectors (a, b, c):\n");
    printf("ax ay az\n");
    printf("bx by bz\n");
    printf("cx cy cz\n");
    printf("\nUsage: %s basisFile gridSpacingX [gridSpacingY gridSpacingZ] outGrid\n", argv[0]);
    return 0;
  }

  const char* basisFile = argv[1];
  double gridSpacingX = strtod(argv[2], NULL);

  // Get system vectors.
  Scatter sysVec(basisFile);
  Matrix3 box = sysVec.topMatrix();
  Grid* orig;
  
  if (argc == 4) orig = new Grid(box, gridSpacingX);
  else {
    double gridSpacingY = strtod(argv[3], NULL);
    double gridSpacingZ = strtod(argv[4], NULL);
    printf("Grid resolution: %g %g %g\n", gridSpacingX, gridSpacingY, gridSpacingZ);
    orig = new Grid(box, gridSpacingX, gridSpacingY, gridSpacingZ);
  }

  char comments[256];
  sprintf(comments, "grid with basis %s ", argv[1]);
  orig->write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  delete orig;
  return 0;
}
