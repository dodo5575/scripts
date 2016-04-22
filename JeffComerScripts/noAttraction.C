#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainRect(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s inGrid outGrid\n", argv[0]);
    return 0;
  }

  Grid g(argv[1]);
  const int n = g.length();

  for (int i = 0; i < n; i++) {
    double v = g.getValue(i);
    if (v < 0.0) g.setValue(i, 0.0);
  }

  char comments[256];
  sprintf(comments, "%s no minima", argv[1]);
  printf("%s\n", comments);
  g.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainRect(argc, argv);
}
