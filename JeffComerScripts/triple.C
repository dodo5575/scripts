#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "useful.H"
#include "Grid.H"

///////////////////////////////////////////////////////////////////////
// Driver

int mainSub(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Usage: %s poreGrid dnaGrid0 dnaGrid1 dnaGrid2 outGrid\n", argv[0]);
    return 0;
  }

  const double dz = 5.0;

  Grid pore(argv[1]);
  Grid dna0(argv[2]);
  Grid dna1(argv[3]);
  Grid dna2(argv[4]);
  
  
  Vector3 o0 = dna0.getOrigin();
  Vector3 o1 = o0 + Vector3(0.0, 0.0, -dz);
  Vector3 o2 = o0 + Vector3(0.0, 0.0, dz);
  Vector3 o3 = o0 + Vector3(0.0, 0.0, -2*dz);
  Vector3 o4 = o0 + Vector3(0.0, 0.0, 2*dz);

  
  dna0.setOrigin(o1);
  dna2.setOrigin(o2);
  

  pore.addInterpolate(dna0);
  pore.addInterpolate(dna1);
  pore.addInterpolate(dna2);
    
  //dna.setOrigin(o3);
  //pore.addInterpolate(dna);
  //dna.setOrigin(o4);
  //pore.addInterpolate(dna);

  char comments[256];
  sprintf(comments, "%s with %s", argv[1], argv[2]);
  printf("%s\n", comments);
  pore.write(argv[argc-1], comments);

  return 0;
}

int main(int argc, char* argv[]) {
  return mainSub(argc, argv);
}
