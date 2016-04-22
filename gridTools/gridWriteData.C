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
  if ( argc != 3 ) {
    printf("Extract the data from the grid and write it as a single column.\n");
    printf("Usage: %s inGrid outDataFile\n", argv[0]);
    return 0;
  }

  const char* inGrid = argv[1];
  const char* outDataFile = argv[argc-1];

  // Load the grid.
  Grid src(inGrid);
  const int size = src.length();
  printf("Loaded `%s' with %d nodes.\n", inGrid, size);


  // Open the file.
  FILE* out = fopen(outDataFile,"w");
  if (out == NULL) {
    printf("ERROR: Couldn't open file `%s'.\n", outDataFile);
    exit(-1);
  }
  
  for (int i = 0; i < size; i++) fprintf(out, "%.12g\n", src.getValue(i));
  fclose(out);
  

  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
