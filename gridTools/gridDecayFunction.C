#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

double decay(double v, double aInf, double a0, double v0) {
  if (v < 0) return a0;
  return aInf*(1.0 - (1.0 - a0/aInf)*exp(-v/v0));
}

int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("v > 0: a(v) = aInf*[1 - (1-a0/aInf)*exp(-v/v0) ]\n");
    printf("v < 0: a(v) = a0\n");
    printf("Usage: %s inGrid aInf a0 v0 outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  double aInf = strtod(argv[2], NULL);
  double a0 = strtod(argv[3], NULL);
  double v0 = strtod(argv[4], NULL);


  const int n = orig.length();
  for (int i = 0; i < n; i++) {
    double a = decay(orig.getValue(i), aInf, a0, v0);
    orig.setValue(i, a);
  }

  char comments[256];
  sprintf(comments, "a(v) = %g*[1 - (1-%g/%g)*exp(-v/%g) ]", aInf, a0, aInf, v0);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
