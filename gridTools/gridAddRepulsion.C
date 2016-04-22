#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
using namespace std;

#include "useful.H"
#include "Grid.H"
#include "CellDecomposition.H"

double energy(double r, double rMinLJ, double vHard) {
  double rHard = 0.45*rMinLJ;
  double rStart = 0.75*rMinLJ;

  if (r >= rStart) return 0.0;
  double dr0 = rHard - rStart;
  double dr = r - rStart;

  return vHard*(dr*dr)/(dr0*dr0);
}

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
  if ( argc != 5 && argc != 6) {
    printf("Add repulsive point charges to a grid.\n");
    printf("Usage: %s inGridFile repulseFile repulseEnergy [radiusScale] outGridFile\n", argv[0]);
    return 0;
  }

  const char* inGridFile = argv[1];
  const char* repulseFile = argv[2];
  const double repulseEnergy = strtod(argv[3], NULL);
  const char* outGridFile = argv[argc-1];

  double radiusScale;
  if (argc == 6) radiusScale = strtod(argv[4], NULL);
  else radiusScale = 1.0;

  // Load the grids.
  Grid g(inGridFile);
  const int size = g.length();
  printf("Loaded grid `%s' of %d nodes.\n", inGridFile, size);
  
  // Load the charges.
  int nSources = countChargeFile(repulseFile);
  double* sourceRad = new double[nSources];
  Vector3* sourcePos = new Vector3[nSources];
  readChargeFile(repulseFile, nSources, sourceRad, sourcePos);
  printf("Read %d points/radii from `%s'.\n", nSources, repulseFile);

  // Scale the radii.
  if (argc == 6)
    for (int i = 0; i < nSources; i++) sourceRad[i] *= radiusScale;

  // Get the largest radius.
  double maxRad = sourceRad[0];
  for (int i = 1; i < nSources; i++) {
    if (sourceRad[i] > maxRad) maxRad = sourceRad[i];
  }

  // Make a cell decomposition using this radius.
  Matrix3 sys = g.getBox();
  Vector3 origin = g.getOrigin();
  printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], maxRad);
  CellDecomposition cell(sys, origin, maxRad);
  printf("Cell decomposition of %i cells\n", cell.length());
  cell.decompose(sourcePos, nSources);
  printf("Cell decomposition contains %i points\n", cell.countPoints());

  // Now go through the grid to create the output potential.
  int complete = -100;
  for (int i = 0; i < size; i++) {
    Vector3 r = g.getPosition(i);
    IndexList neigh = cell.neighborhood(r);
   
    // Go through the neighboring sources.
    int nNeighs = neigh.length();
    double sum = 0.0;
    for (int j = 0; j < nNeighs; j++) {
      int ind = neigh.get(j);
      Vector3 d = g.wrapDiffNearest(sourcePos[ind]-r);
      double dist = d.length();
      sum += energy(dist, sourceRad[ind], repulseEnergy);
    }

    // Add the energy to the grid point.
    g.setValue(i, g.getValue(i)+sum);

    // Write the progress.
    int comp = (100*i)/size;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  }
  printf("Added the potential energy due to the sources.\n");
  
    // Write the result.
  char comments[256];
  sprintf(comments, "%s with %s sources, repulseEnergy=%g", inGridFile, repulseFile, repulseEnergy);
  g.write(outGridFile, comments);
  printf("%s\n", comments);
  printf("Wrote `%s'.\n", outGridFile);
  
  return 0;
}

int main(int argc, char* argv[]) {
  return mainAdd(argc, argv);
}
