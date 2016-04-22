#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Blur the data blurCount times.\n");
    printf("Usage: %s inGrid maskGrid blurCount outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  Grid mask(argv[2]);
  int blurCount = atoi(argv[3]);

  // Do some blurring.
  for (int i = 0; i < blurCount; i++) orig.blur(mask);

  char comments[256];
  snprintf(comments, 256, "%s blurred %d times", argv[1], blurCount);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
