#include <cmath>
#include <cstdio>

using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int main(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Add a cylindrical pore with a barrier potential of `potential'.\n");
    printf("Usage: %s inGrid poreDiameter poreLength potential outGrid\n", argv[0]);
    return 0;
  }

  Grid orig(argv[1]);
  double poreDiameter = strtod(argv[2], NULL);
  double poreLength = strtod(argv[3], NULL);
  double potential = strtod(argv[3], NULL);

  double ss0 = 0.25*(poreDiameter*poreDiameter);
  double l0 = 0.5*poreLength;

  const int n = orig.getSize();
  for (int i = 0; i < n; i++) {
    Vector3 r = orig.getPosition(i);
    
    if (fabs(r.z) < l0) {
      double ss = r.x*r.x + r.y*r.y;
      double v = orig.getValue(i);

      if (ss > ss0) orig.setValue(i, potential+v);
    }
  }

  char comments[256];
  sprintf(comments, "cylindrical pore, d=%g, l=%g", poreDiameter, poreLength);
  orig.write(argv[argc-1], comments);
  printf("Wrote `%s'.\n", argv[argc-1]);

  return 0;
}
