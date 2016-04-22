#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Usage: %s srcFile z0 z1 destFile\n", argv[0]);
    return 0;
  }

  double z0 = strtod(argv[2], NULL);
  double z1 = strtod(argv[3], NULL);

  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  const int n = src.length();
  
  int count = 0;
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);
    if (r.z >= z0 && r.z < z1) {
      sum += src.getValue(i);
      count++;
    }
  }
  
  double avg = sum/count;
  printf("Average between z = %.4g and z = %.4g:\n", z0, z1);
  printf("%.10g\n", avg);

  // Shift the grid to make the average zero.
  src.shift(-avg);
  
  char comments[256];
  sprintf(comments, "shifted %s by %.4g", argv[1], -avg);
  printf("%s\n", comments);
  src.write(argv[argc-1], comments);

  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
