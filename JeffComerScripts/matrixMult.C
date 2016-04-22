#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"

///////////////////////////////////////////////////////////////////////
// Driver


int mainAdd(int argc, char* argv[]) {
  if ( argc < 19 ) {
    printf("Usage: %s axx axy axz ayx ayy ayz azx azy azz bxx bxy...\n", argv[0]);
    return 0;
  }
  double a[9];
  double b[9];


  for (int i = 0; i < 9; i++) {
    a[i] = atoi(argv[i+1]);
    b[i] = atoi(argv[i+10]);
  }
  
  Matrix3 ma(a);
  Matrix3 mb(b);
  
  Matrix3 mc = ma*mb;
  printf("%s\n", mc.toString().val());

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
