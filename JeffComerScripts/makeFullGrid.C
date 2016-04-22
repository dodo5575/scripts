#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver
int mainFullGrid(int argc, char* argv[]) {
  if ( argc != 7 ) {
    printf("Usage: %s partGrid lx ly lz dl outGrid\n", argv[0]);
    return 0;
  }

  Grid part(argv[1]);
  printf("Read %s containing %d nodes.\n", argv[1], part.length());
  Vector3 l;
  l.x = strtod(argv[2],NULL);
  l.y = strtod(argv[3],NULL);
  l.z = strtod(argv[4],NULL);
  double dl = strtod(argv[5],NULL);

  // Get the dimensions of part.
  Vector3 p0 = part.getOrigin();
  Vector3 p1 = part.getDestination();


  printf("Generating new grid with dimensions of %s.\n", l.toString().val());
  Grid full(l, dl);
  Vector3 newOrigin(0.0, 0.0, 0.5*dl);
  full.setOrigin(full.getOrigin() + newOrigin);

  const int n = full.length();
  for (int i = 0; i < n; i++) {
    Vector3 r = full.getPosition(i);

    // Reflect the PMF about z = 0.
    if (r.z < p0.z + 3*dl) r.z = -r.z;

    // Clone the sides.
    if (r.x < p0.x) r.x = p0.x;
    if (r.x >= p1.x) r.x = p1.x;
    if (r.y < p0.y) r.y = p0.y;
    if (r.y >= p1.y) r.y = p1.y;
    if (r.z < p0.z + 3*dl) r.z = p0.z - 3*dl;
    if (r.z >= p1.z - 3*dl) r.z = p1.z - 3*dl;
    
    // Interpolate the current position.
    full.setValue(i, part.interpolatePotential(r));
    //full.setValue(i, part.getPotential(r));
    //printf("pos: %s\n", r.toString().val());
  }

  char comments[256];
  sprintf(comments, "Extended %s to fill %s.", argv[1], l.toString().val());
  full.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainFullGrid(argc, argv);
}
