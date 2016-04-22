#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainProfileReflect(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s inGrid externalForce\n", argv[0]);
    return 0;
  }

  double externalForce = strtod(argv[2], NULL);

  char outNeg[256];
  char outPos[256];
  sprintf(outNeg, "%s.neg.dat", argv[1]);
  sprintf(outPos, "%s.pos.dat", argv[1]);
  
  Grid pos(argv[1]);
  Grid neg(pos);

  Vector3 force(0.0, 0.0, -externalForce);
  double v0;

  pos.addGradient(force);
  v0 = pos.getValue(0);
  pos.shift(-v0);
    

  neg.flipZ();
  neg.addGradient(force);
  v0 = pos.getValue(0);
  neg.shift(-v0);
  
  pos.averageProfileZBoltzmann(outPos);
  neg.averageProfileZBoltzmann(outNeg);
}

int main(int argc, char* argv[]) {
  return mainProfileReflect(argc, argv);
}
