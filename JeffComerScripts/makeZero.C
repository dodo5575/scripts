#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainGate(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s inGrid outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  Grid out(orig);
  out.zero();
  char comments[256];
  sprintf(comments, "zeroed %s",argv[1]);
  printf("%s\n", comments);
  out.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainGate(argc, argv);
}
