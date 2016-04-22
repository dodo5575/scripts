#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  char usage[] = "Usage: %s [-i] inGrid threshold outGrid\n";

  if ( argc != 4 && argc != 5) {
    printf("All values < threshold are set to 0. All values >= threshold are set to 1.\n");
     printf("-i signifies that the result is the inverse mask: v'(i) = 1-v(i)\n");
    printf(usage, argv[0]);
    return 0;
  }
  
  // Check for inversion flag.
  bool invert = false;
  if (argc == 5) {
    if (argv[1][0] == '-' && argv[1][1] == 'i') invert = true;
    else {
      printf(usage, argv[0]);
      return 0;
    }
  }

  const char* inGrid = argv[argc-3];
  double threshold = strtod(argv[argc-2], NULL);
  const char* outGrid = argv[argc-1];

  Grid orig(inGrid);
  orig.alphaThreshold(threshold);
  
  // Invert if necessary.
  if (invert) orig.alphaInvert();

  char comments[256];
  snprintf(comments, 256, "%s alpha threshold %g", inGrid, threshold);
  orig.write(outGrid, comments);
  printf("Wrote `%s'.\n", outGrid);

  return 0;
}
