#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("All values < minValue are set to minValue. All values > maxValue are set to maxValue.\n");
    printf("Usage: %s inGrid minValue maxValue outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  double minValue = strtod(argv[2], NULL);
  double maxValue = strtod(argv[3], NULL);

  orig.threshold(minValue, maxValue);

  char comments[256];
  sprintf(comments, "%s threshold %g %g", argv[1], minValue, maxValue);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
