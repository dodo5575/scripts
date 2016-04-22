#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainDo(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s gridSize outGrid\n", argv[0]);
    return 0;
  }

  double gridSize = strtod(argv[1], NULL);

  Vector3 lx = Vector3(66.7379, 38.539, 0);
  Vector3 ly = Vector3(0, 77.0781, 0);
  Vector3 lz = Vector3(0, 0, 227.994);
  Matrix3 box(lx, ly, lz);

  
  Grid g(box, gridSize);
  char comment[256];
  sprintf(comment, "grid size: %g", gridSize); 
  g.write(argv[2], comment);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainDo(argc, argv);
}
