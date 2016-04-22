#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
//#include "Scatter.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Usage: %s srcFile z0 blurCount outputFile\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  
  double z0 = strtod(argv[2],NULL);
  int blurCount = atoi(argv[3]);

  Vector3 ez = src.getBasis().ez();
  Vector3 origin = src.getOrigin();
  //src.setOrigin(origin + Vector3(0.0, 0.0, 0.5*ez.z));

  // Zero everything within farther along the z axis than z0.
  const int n = src.length();
  for (int i = 0; i < n; i++) {
    Vector3 pos = src.getPosition(i);
    if (fabs(pos.z) > z0) src.setValue(i,0.0);
  }

  // Blur the specified number of times.
  for (int i = 0; i < blurCount; i++) src.blur();

  src.write(argv[argc-1], "touched up");
  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
