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
  if ( argc < 4 ) {
    printf("Generate a boundary grid from a mask.\n");
    printf("Usage: %s maskGrid blurCount outGrid\n", argv[0]);
    return 0;
  }

  const int blurCount = atoi(argv[2]);

  // Make a boundary mask.
  Grid m(argv[1]);
  Grid bord(m);
  bord.alphaThreshold(0.6);
  for (int i = 0; i < blurCount; i++) bord.blur();
  Grid inv(bord);
  inv.alphaInvert();
  bord.multiply(inv);
  bord.alphaThreshold(0.0001);
  bord.multiply(m);
  bord.alphaThreshold(0.0001);
  bord.blur();

  // Write the result.
  char comments[256];
  sprintf(comments, "%s border", argv[1]);
  printf("%s\n", comments);
  bord.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
