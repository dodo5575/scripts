// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"
#include "Scatter.H"

using namespace std;

int countCoordRadius(const char* fileName) {
  int nRead;
  int n = 0;
  double x, y, z, rad;
  char line[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("countCoordRadius Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%lf %lf %lf %lf", &rad, &x, &y, &z);
    if (nRead >= 4) n++;
  }
    
  fclose(inp);
  return n;
}

int readCoordRadius(const char* fileName, Vector3* pos, double* radius) {
  int nRead;
  int n = 0;
  double x, y, z, rad;
  char line[256];

  // Open the file.
  FILE* inp = fopen(fileName,"r");
  if (inp == NULL) {
    printf("readCoordRadius Couldn't open file %s\n.",fileName);
    exit(-1);
  }

  while (fgets(line, 256, inp) != NULL) {
    // Ignore comments.
    int len = strlen(line);
    if (line[0] == '#') continue;
    if (len < 2) continue;
    
    // Read values.
    nRead = sscanf(line, "%lf %lf %lf %lf", &rad, &x, &y, &z);
    if (nRead >= 4) {
      pos[n] = Vector3(x, y, z);
      radius[n] = rad;
      n++;
    } else {
      printf("readCoordRadius needs 'radius x y z' format.\n");
      exit(-1);
    }
  }
    
  fclose(inp);
  return n;
}


///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
  if ( argc != 5 ) {
    printf("Set the grid to distances from the point.\n");
    printf("Usage: %s srcGrid radiusPositionFile srcRadius outFile\n", argv[0]);
    return 0;
  }

  const char* srcFile = argv[1];
  const char* pointFile = argv[2];
  const double srcRadius = strtod(argv[3], NULL);
  const char* outFile = argv[argc-1];

  // Load the grids.
  Grid src(srcFile);
  src.zero();
  int n = src.length();
  printf("Loaded `%s' which contains %d nodes.\n", srcFile, n); 

  // Load the coordinates.
  printf("Loading the coordinates.\n");
  int posNum = countCoordRadius(pointFile);
  Vector3* pos = new Vector3[posNum];
  double* radius = new double[posNum];
  
  readCoordRadius(pointFile, pos, radius);
  printf("Loaded %d points from `%s'.\n", posNum, pointFile);

  int complete = 0;
  for (int i = 0; i < n; i++) {
    Vector3 r = src.getPosition(i);

    Vector3 dr = pos[0]-r;
    double minDist = dr.length() - radius[0] - srcRadius;
    for (int j = 1; j < posNum; j++) {
      Vector3 dr = pos[j] - r;
      double dist = dr.length() - radius[j] - srcRadius;

      if (dist < minDist) minDist = dist;
    }
    src.setValue(i, minDist);

    // Write the progress.
    int comp = (100*i)/n;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  }
  char comments[256];
  sprintf(comments, "%s patched", srcFile);
  src.write(outFile, comments);
  printf("Wrote `%s'.\n", outFile);

  delete[] pos;
  delete[] radius;

  return 0;
}
