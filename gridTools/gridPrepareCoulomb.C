#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"


int countChargeFile(const char* fileName) {
  int nRead;
  int n = 0;
  double q, x, y, z;
  char line[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%lf %lf %lf %lf", &q, &x, &y, &z);
    if (nRead >= 4) n++;
  }
    
  fclose(inp);
  return n;
}

// Read coordinates into a Vector array.
void readChargeFile(const char* fileName, int num, double* charge, Vector3* pos) {
  int nRead;
  int n = 0;
  double q, x, y, z;
  char line[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("Scatter:countCoordinates Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%lf %lf %lf %lf", &q, &x, &y, &z);
    if (nRead >= 4) {
      charge[n] = q;
      pos[n].x = x;
      pos[n].y = y;
      pos[n].z = z;
      n++;
    }
  }
    
  fclose(inp);
}
///////////////////////////////////////////////////////////////////////
// Driver

int mainAdd(int argc, char* argv[]) {
  if ( argc != 6 ) {
    printf("Add one grid to another applying the transformation given in gridListFile.\n");
    printf("Usage: %s inGridFile maskGridFile coulombFile coulombConstant outGridFile\n", argv[0]);
    printf("coulombFile format: 'q x y z'\n");
    return 0;
  }

  const char* inGridFile = argv[1];
  const char* maskGridFile = argv[2];
  const char* coulombFile = argv[3];
  const double coulombConst = strtod(argv[4], NULL);
  const char* outGridFile = argv[argc-1];

  char bordGridFile[256];
  sprintf(bordGridFile, "%s.bord", inGridFile);

  const int blurLevel = 10;

  // Load the grids.
  Grid src(inGridFile);
  printf("Loaded grid of %d nodes.\n", src.length());
  Grid mask(maskGridFile);
  printf("Loaded mask of %d nodes.\n", mask.length());
  // Make sure that this is actually a mask.
  mask.threshold(0.0, 1.0);

  // Load the charges.
  int nCharges = countChargeFile(coulombFile);
  double* charge = new double[nCharges];
  Vector3* chargePos = new Vector3[nCharges];
  readChargeFile(coulombFile, nCharges, charge, chargePos);
  printf("Read %d charges from %s.\n", nCharges, coulombFile);
  
  // Subtract the charges.
  for (int p = 0; p < nCharges; p++) src.addCoulomb(-coulombConst*charge[p], chargePos[p]);
  printf("Subtracted the potential energy due to the charges.\n");
  
  // Make a boundary mask.
  Grid bord(mask);
  bord.alphaThreshold(0.6);
  for (int i = 0; i < blurLevel; i++) bord.blur();
  Grid inv(bord);
  inv.alphaInvert();
  bord.multiplyInterpNoWrap(inv);
  bord.alphaThreshold(0.01);
  bord.multiplyInterpNoWrap(mask);
  bord.alphaThreshold(0.99);
  bord.blur();
  bord.write(bordGridFile);

  // Zero based on the average of the boundary nodes.
  double v0 = src.averageRegion(bord);
  printf("Average potential energy on border: %g\n", v0);
  src.shift(-v0);

  // Apply the mask.
  src.multiplyInterpNoWrap(mask);
 
  // Write the result.
  char comments[256];
  sprintf(comments, "%s sum", inGridFile);
  src.write(outGridFile, comments);
  printf("%s\n", comments);
  printf("Wrote `%s'.\n", outGridFile);
  
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
