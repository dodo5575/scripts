// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Create a grid = exp(-grid)\n");
    printf("Usage: %s srcGrid outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const char* outFile = argv[argc-1];

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded `%s'.\n", inFile); 
  src.boltzmann();

  char comments[256];
  sprintf(comments, "%s Boltzmann", inFile);
  src.write(outFile, comments);

  return 0;
}
