#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 10 ) {
    printf("Usage: %s inGrid Xo Yo Zo R1 T1 R2 T2 outGrid\n", argv[0]);
    printf("\nGenerates temperature map for hot sphere [T = T1 @ R = R1] with boundary conditions [T = T2 @ R = R2].\n");
    printf("r > R1: T(r) = A/r + B; where sphere is centered at r0 = (Xo, Yo, Zo).\n");
    printf("        A = (T1 - T2) * R1 * R2 / (R2 - R1)\n");
    printf("        B = (T2 * R2 - T1 * R1) / (R2 - R1)\n");
    printf("r < R1: T(r) = T(R1)\n");

    return 0;
  }

  const char* inFile  = argv[1];
  const char* outFile = argv[argc-1];
  double Xo = strtod(argv[2], NULL);
  double Yo = strtod(argv[3], NULL);
  double Zo = strtod(argv[4], NULL);
  double R1 = strtod(argv[5], NULL);
  double T1 = strtod(argv[6], NULL);
  double R2 = strtod(argv[7], NULL);
  double T2 = strtod(argv[8], NULL);

  Grid src(inFile);
  printf("Loaded `%s'.\n", inFile);
  src.tempSphere(Vector3(Xo,Yo,Zo), R1, T1, R2, T2);

  char comments[256];
  sprintf(comments, "Hot Sphere: R1=%g, T1=%g, R2=%g, T2=%g",R1,T1,R2,T2);
  src.write(outFile, comments); 

  return 0;
}
