#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainSub(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s inGrid outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  Grid att(orig);
  const int n = orig.length();

  for (int i = 0; i < n; i++) {
    double v = att.getValue(i);
    if (v > -0.3) att.setValue(i, 0.0);
  }

  char comments[256];
  sprintf(comments, "%s remove", argv[1]);
  printf("%s\n", comments);
  att.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainSub(argc, argv);
}
