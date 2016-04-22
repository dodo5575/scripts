#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 4 ) {
    printf("Blur the data blurCount times.\n");
    printf("Usage: %s inGrid blurCount outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  int blurCount = atoi(argv[2]);

  // Do some blurring.
  for (int i = 0; i < blurCount; i++) orig.blur();

  char comments[256];
  sprintf(comments, "%s blurred %d times", argv[1], blurCount);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
