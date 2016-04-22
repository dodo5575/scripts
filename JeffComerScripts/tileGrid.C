#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Usage: %s srcFile nx ny nz outputFile\n", argv[0]);
    return 0;
  }

  // Load the first grid.
  Grid src(argv[1]);
  printf("Loaded %s.\n", argv[1]);

  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[4]);
  printf("Tiling %d by %d by %d.\n", nx, ny, nz);

  Grid dest(src.tile(nx, ny, nz));

  char comments[256];
  sprintf(comments, "%s tiled %d by %d by %d", argv[1], nx, ny, nz);
  dest.write(argv[argc-1], comments);
  printf("%s\n", comments);

  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
