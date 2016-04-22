#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Homogenize along z.\n");
    printf("Usage: %s inGrid outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  Grid dest(orig.threeDToOneD());

  char comments[256];
  snprintf(comments, 256, "%s homogenized", argv[1]);
  dest.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
