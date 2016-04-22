#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainCrossSection(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Generate 2D cross sections of the grid data.\n");
    printf("Usage: %s inGrid dir factor outData\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  int dir = atoi(argv[2]);
  double factor = strtod(argv[3],NULL);

  // Do some blurring.
  //for (int i = 0; i < 3; i++) orig.blur();

  orig.crossSectionFactor(argv[argc-1], dir, factor);
  printf("Wrote %s.\n", argv[argc-1]);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainCrossSection(argc, argv);
}
