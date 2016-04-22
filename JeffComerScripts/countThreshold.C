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
  if ( argc != 3 ) {
    printf("Usage: %s srcFile threshold\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  
  double threshold = strtod(argv[2],NULL);

  int count = 0;
  const int n = src.length();
  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    if (v < threshold) count++;
  }

  printf("total: %i\n", n);
  printf("count: %i\n", count);
  printf("fraction: %.10g\n", double(count)/n);
  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
