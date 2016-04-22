#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainGate(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Usage: %s inGrid z0 len height outGrid\n", argv[0]);
    return 0;
  }
  const double pi = 4.0*atan(1.0);

  Grid orig(argv[1]);
  const double z0 = strtod(argv[2], NULL);
  const double len = strtod(argv[3], NULL);
  const double height = strtod(argv[4], NULL);

  const double z1 = -z0 - len;
  const double z2 = -z0;
  const double z3 = z0;
  const double z4 = z0 + len;

  const int n = orig.length();
  for (int i = 0; i < n; i++) {
    Vector3 r = orig.getPosition(i);

    if (fabs(r.z) > z4) continue;
    double v0 = orig.getValue(i);
    double v;
    
    if (r.z >= z1 && r.z < z2) {
      v = 0.5*height*(1 + cos(pi*(r.z+z0)/len));
    } else if (r.z >= z2 && r.z < z3) {
      v = height;
    } else {
      v =  0.5*height*(1 + cos(pi*(r.z-z0)/len));
    }
    orig.setValue(i, v + v0);
  }

  char comments[256];
  sprintf(comments, "%s gated", argv[1]);
  printf("%s\n", comments);
  orig.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainGate(argc, argv);
}
