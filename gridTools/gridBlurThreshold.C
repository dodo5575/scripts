#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Blur the data blurCount times for points with values bigger than threshold.\n");
    printf("Usage: %s inGrid blurCount threshold outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  int blurCount = atoi(argv[2]);
  double threshold = strtod(argv[3], NULL);

  // Do some blurring.
  for (int i = 0; i < blurCount; i++) orig.blurBigger(threshold);

  char comments[256];
  sprintf(comments, "%s blurred %d times", argv[1], blurCount);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
