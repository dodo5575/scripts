#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s srcFile destFile\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  const int n = src.length();
  
  // Get the minimum value of the density that is not zero.
  double minVal = src.getValue(0);
  for (int i = 1; i < n; i++) {
    double v = src.getValue(i);
    if (v > 0.0 && v < minVal) minVal = v;
  }
  
  // -log(minVal) will stand in for infinity.
  double infinity = -log(minVal);
  
  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    if (v <= 0.0) src.setValue(i, infinity);
    else src.setValue(i, -log(v));
  }

  char comments[256];
  sprintf(comments, "-log of %s", argv[1]);
  printf("%s\n", comments);
  src.write(argv[argc-1], comments);

  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
