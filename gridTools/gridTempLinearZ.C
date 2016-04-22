#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 9 ) {
    printf("Usage: %s inGrid z1 z2 t1 z3 z4 t2 outGrid\n", argv[0]);
    printf("Generates temperature map for the system where:.\n");
    printf("    T(z1 < z < z2) = T1\n");
    printf("    T(z3 < z < z4) = T2\n");
    printf("And a linear fit is used between the regions.\n");

    return 0;
  }

  const char* inFile  = argv[1];
  const char* outFile = argv[argc-1];
  double z1 = strtod(argv[2], NULL);
  double z2 = strtod(argv[3], NULL);
  double t1 = strtod(argv[4], NULL);
  double z3 = strtod(argv[5], NULL);
  double z4 = strtod(argv[6], NULL);
  double t2 = strtod(argv[7], NULL);

  Grid src(inFile);
  printf("Loaded `%s'.\n", inFile);
  src.tempGradZ(z1,z2,t1,z3,z4,t2);

  char comments[256];
  sprintf(comments, "System with two temperature-controlled regions: T(%g<z<%g)=%g, T(%g<z<%g)=%g",z1,z2,t1,z3,z4,t2);
  src.write(outFile, comments); 

  return 0;
}
