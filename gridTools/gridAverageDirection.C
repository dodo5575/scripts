// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

double min(double x, double y) {
  if (x <= y) return x;
  return y;
}

double max(double x, double y) {
  if (x >= y) return x;
  return y;
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Average the grid along a given axis.\n");
    printf("Usage: %s srcGrid direction(x|y|z) outFile\n", argv[0]);
    return 0;
  }

  const char* inFile = argv[1];
  const char dir = argv[2][0];
  const char* outFile = argv[argc-1];

  int direct = -1;
  if (dir == 'X' || dir == 'x') direct = 0;
  if (dir == 'Y' || dir == 'y') direct = 1;
  if (dir == 'Z' || dir == 'z') direct = 2;
  
  if (direct < 0) {
    printf("Error: Direction must be x, y, or z!\n");
    return 0;
  }

  // Load the first grid.
  Grid src(inFile);
  printf("Loaded `%s'.\n", inFile); 
  printf("direction: %d\n", direct);
  
  src.average(direct);
  char comments[256];
  snprintf(comments, 256, "%s averaged along %c", argv[1], dir);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);
  
  return 0;
}
