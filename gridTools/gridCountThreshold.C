#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "useful.H"
//#include "Scatter.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int mainTouch(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Count number of nodes less than the threshold.\n");
    printf("Usage: %s srcFile threshold\n", argv[0]);
    return 0;
  }

  const double waterVol = 29.91; // angstrom^3

  // Load the first grid.
  Grid src(argv[1]);
  printf("\nLoaded %s.\n", argv[1]);
  
  double threshold = strtod(argv[2],NULL);

  int count = 0;
  const int n = src.length();
  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    if (v < threshold) count++;
  }
  double vol = src.getVolume();
  double frac = double(count)/n;

  printf("threshold: %.10g\n", threshold);
  printf("total: %i\n", n);
  printf("count: %i\n", count);
  printf("fraction: %.10g\n", frac);
  printf("\nvolume_total: %.10g\n", vol);
  printf("volume: %.10g\n", frac*vol);
  printf("water_count: %.10g\n", frac*vol/waterVol);

  return 0;
}


int main(int argc, char* argv[]) {
  return mainTouch(argc, argv);
}
