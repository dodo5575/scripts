#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

void printFormat() {
  printf("transformFile contains the transformation matrix and optionally a displacement:\n");
  printf("exx exy exz\n");
  printf("eyx eyy eyz\n");
  printf("ezx ezy ezz\n");
  printf("[dx dy dz]\n");
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 4 && argc != 5 ) {
    printf("Transform the grid using a matrix and a displacement.\n");
    printFormat();
    printf("\nUsage: %s [-i] inGrid tranformFile outGrid\n", argv[0]);
    printf("-i inverts the transformation.\n");
    return 0;
  }
  
  const char* inGrid = argv[argc-3];
  const char* transFile = argv[argc-2];
  const char* outGrid = argv[argc-1];

  // Invert the transformation?
  bool invert = false;
  if (argc == 5 && argv[1][0] == '-' && argv[1][1] == 'i') invert = true;

  // Get system vectors.
  Scatter trans(transFile);
  if (trans.length() < 3) {
    printf("gridTransform Transform file `%s' does not contain enough entries.\n", transFile);
    printFormat();
    exit(-1);
  }
  if (trans.length() > 4) {
    printf("gridTransform Transform file `%s' contains too many entries.\n", transFile);
    printFormat();
    exit(-1);
  }
  
  // Get the transformation matrix.
  Matrix3 t = trans.rowMatrix();
  Matrix3 it = t.inverse();
  // Get the displacement if it exists.
  Vector3 disp(0.0);
  if (trans.length() == 4) disp = trans.get(3);

  // Load the grid.
  Grid src(inGrid);
  printf("Loaded `%s'.\n", inGrid);

  // Transform it.
  Vector3 origin = src.getOrigin();
  Matrix3 basis = src.getBasis();
  if (invert) {
    origin = it.transform(origin + disp);
    basis = it.transform(basis);
  } else {
    origin = t.transform(origin) + disp;
    basis = t.transform(basis);
  }

  // Set the results.
  src.setOrigin(origin);
  src.setBasis(basis);

  char comments[256];
  snprintf(comments, 256, "grid %s transformed by %s", inGrid, transFile);
  src.write(outGrid, comments);
  printf("%s\n", comments);
  printf("Wrote `%s'.\n", outGrid);

  return 0;
}
