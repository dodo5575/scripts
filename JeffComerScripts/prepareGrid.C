#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s srcFile templateFile outputFile\n", argv[0]);
    return 0;
  }
  const double minValue = 20.0;

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded %s.\n", argv[1]);

  Grid temp(argv[2]);
  printf("Loaded %s.\n", argv[2]);

  printf("%s\n", temp.getBox().toString().val());
  Grid dest(temp.getBox(), src.getNx(), src.getNy(), src.getNz());

  const int n = dest.getSize();
  for (int i = 0; i < n; i++) {
    Vector3 r = dest.getPosition(i);
    double v = src.interpolatePotential(r);
    if (v < minValue) dest.setValue(i, minValue);
    else dest.setValue(i, v);
  }

  char comments[256];
  sprintf(comments, "prepared diffusion from %s", argv[1]);
  dest.write(argv[argc-1], comments);
  printf("%s\n", comments);

  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
