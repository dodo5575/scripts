#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainSub(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s inGrid0 inGrid1 outGrid\n", argv[0]);
    return 0;
  }

  Grid g0(argv[1]);
  Grid g1(argv[2]);
  
  g0.add(g1);
  g0.scale(0.5);

  char comments[256];
  sprintf(comments, "%s - %s", argv[1], argv[2]);
  printf("%s\n", comments);
  g0.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainSub(argc, argv);
}
