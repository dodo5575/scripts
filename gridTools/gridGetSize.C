// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"

using namespace std;

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 2 ) {
    printf("Usage: %s srcGrid\n", argv[0]);
    return 0;
  }

  const char* srcGrid = argv[1];
  
  Grid src(srcGrid);
  Vector3 origin = src.getOrigin();
  Vector3 dest = src.getDestination();
  Vector3 center = 0.5*(origin+dest);
  Vector3 ext = src.getExtent();
  Vector3 ex = src.getBasis().ex();
  Vector3 ey = src.getBasis().ey();
  Vector3 ez = src.getBasis().ez();
  Matrix3 box = src.getBox();

  printf("%s\n", src.getBasis().toString().val());
  printf("nodes %d %d %d\n", src.getNx(), src.getNy(), src.getNz());
  printf("origin %.10g %.10g %.10g\n", origin.x, origin.y, origin.z);
  printf("\n");
  printf("size %d\n", src.length());
  printf("extent %.10g %0.10g %.10g\n", ext.x, ext.y, ext.z);
  printf("diagonal %.10g\n", ext.length());
  printf("volume %.10g\n", src.getVolume());
  printf("a %.10g\n", ex.length());
  printf("b %.10g\n", ey.length());
  printf("c %.10g\n", ez.length());
  printf("center %.10g %.10g %.10g\n", center.x, center.y, center.z);
  printf("\n");  
  printf("cellBasisVector1 %s\n", box.ex().toString().val());
  printf("cellBasisVector2 %s\n", box.ey().toString().val());
  printf("cellBasisVector3 %s\n", box.ez().toString().val());

  return 0;
}
