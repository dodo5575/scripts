#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainAdd(int argc, char* argv[]) {
  if ( argc < 3 ) {
    printf("Usage: %s outGrid inGrid0 inGrid1...\n", argv[0]);
    return 0;
  }

  const int n = argc-2;
  printf("Adding %d grids.\n", n);
  
  // Load the first grid.
  Grid g0(argv[2]);

  // Interpolate the others onto it.
  for (int i = 3; i < argc; i++) {
    Grid g(argv[i]);
    g0.addLinear(g);
  }

  // Write the result.
  char comments[256];
  sprintf(comments, "%s sum", argv[2]);
  printf("%s\n", comments);
  g0.write(argv[1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
