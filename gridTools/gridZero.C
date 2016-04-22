#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Set all grid entries to zero.\n");
    printf("Usage: %s inGrid outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  orig.zero();

  char comments[256];
  sprintf(comments, "%s zeroed", argv[1]);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
