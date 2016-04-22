#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
// 
// Makes an emtpy DX file.  Great for defining periodic boundaries
// and bins for the 3D WHAM PMF calculations.
// 
///////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
  char s[128];

  if ( argc != 11 ) {
    printf("Usage: %s outFile ox oy oz lx ly lz nx ny nz\n(you gave %d arguments)\n", argv[0],argc);
    return 0;
  }
  
  Vector3 origin(strtod(argv[2], NULL),strtod(argv[3], NULL),strtod(argv[4], NULL));
  Matrix3 box(strtod(argv[5], NULL)/strtod(argv[8], NULL),strtod(argv[6], NULL)/strtod(argv[9], NULL),strtod(argv[7], NULL)/strtod(argv[10], NULL));
  Grid g(box, origin, atoi(argv[8]), atoi(argv[9]), atoi(argv[10]) );
  

  g.write(argv[1], "");

  printf("\nWrote %s.\n", argv[1]);

  return 0;

 }
