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
  if ( argc != 5 ) {
    printf("Add Coulomb point charges to a grid.\n");
    printf("Usage: %s inGridFile coulombFile coulombConstant outGridFile\n", argv[0]);
    return 0;
  }

  const char* inGridFile = argv[1];
  const char* coulombFile = argv[2];
  const double coulombConst = strtod(argv[3], NULL);
  const char* outGridFile = argv[argc-1];

  // Load the grids.
  Grid src(inGridFile);
  printf("Loaded grid of %d nodes.\n", src.length());
  
  // Load the charges.
  int nCharges = countChargeFile(coulombFile);
  double* charge = new double[nCharges];
  Vector3* chargePos = new Vector3[nCharges];
  readChargeFile(coulombFile, nCharges, charge, chargePos);
  printf("Read %d charges from %s.\n", nCharges, coulombFile);
  
  // Add the charges.
  int complete = -100;
  for (int p = 0; p < nCharges; p++) {
    src.addCoulomb(coulombConst*charge[p], chargePos[p]); 

    // Write the progress.
    int comp = (100*p)/nCharges;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  }
  printf("Added the potential energy due to the charges.\n");
  
    // Write the result.
  char comments[256];
  sprintf(comments, "%s with coulomb charges, coulombConst=%g", inGridFile, coulombConst);
  src.write(outGridFile, comments);
  printf("%s\n", comments);
  printf("Wrote `%s'.\n", outGridFile);
  
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
