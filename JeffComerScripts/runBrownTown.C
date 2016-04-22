///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>    
#include <stdlib.h>
#include "GrandBrownTown.H"

int main(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Usage: %s configFile outputName randomSeed\n", argv[0]);
    printf("You entered %i arguments.\n", argc-1);
    return 0;
  }

  const char* configFile = argv[1];
  const char* outArg = argv[2];
  const long int randomSeed = atol(argv[3]);

  GrandBrownTown brown(configFile, outArg, randomSeed);
  printf("Running...\n");
  brown.run();

  return 0;
}
