#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainRect(int argc, char* argv[]) {
  if ( argc != 3 ) {
    printf("Usage: %s inGrid outGrid\n", argv[0]);
    return 0;
  }

  Grid g(argv[1]);
  
  g.uniformSection(Vector3(0.0));

  const int n = g.length();
  double radius0 = 2;
  double rSq = radius0*radius0;

  for (int i = 0; i < n; i++) {
    Vector3 r = g.getPosition(i);
    if (r.x*r.x + r.y*r.y < rSq) g.setValue(i, 0.0);
  }

  char comments[256];
  sprintf(comments, "%s as a tube", argv[1]);
  printf("%s\n", comments);
  g.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainRect(argc, argv);
}
