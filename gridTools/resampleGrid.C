#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainAdd(int argc, char* argv[]) {
  if ( argc < 4 ) {
    printf("Usage: %s inGrid nx ny nz blurCount outGrid\n", argv[0]);
    return 0;
  }
  
  Grid grid(argv[1]);
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[4]);
  int blurCount = atoi(argv[5]);
  
  grid.resample(nx, ny, nz);
  for (int i = 0; i < blurCount; i++) grid.blur();

  char comments[256];
  sprintf(comments, "%s resampled", argv[1]);
  grid.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
