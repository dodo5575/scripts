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
    printf("Usage: %s srcFile outputFile\n", argv[0]);
    return 0;
  }
  
  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);

  src.reflectZ();
  src.write(argv[argc-1], "mirrored");

  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
