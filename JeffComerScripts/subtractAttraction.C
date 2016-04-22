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
    if (v > 0.0) att.setValue(i, 0.0);
    else att.setValue(i, -att.getValue(i));
  }

  char comments[256];
  sprintf(comments, "%s remove", argv[1]);
  printf("%s\n", comments);
  att.write(argv[argc-1], comments);

  return 0;
}

int mainAdd(int argc, char* argv[]) {
  if ( argc < 2 ) {
    printf("Usage: %s addGrid grid1 grid2...\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  for (int i = 2; i < argc; i++) {
    Grid sum(orig);
    Grid g(argv[i]);
    sum.add(g);
    
    char fileName[256];
    sprintf(fileName, "harder_%s", argv[i]);
    char comments[256];
    sprintf(comments, "%s hard", argv[1]);

    sum.write(fileName, comments);
    printf("%s\n", fileName);
  }

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
