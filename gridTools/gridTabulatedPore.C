// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <omp.h>

#include "useful.H"
#include "Grid.H"
#include "CylinderPore.H"
#include "CellDecomposition.H"
#include "TabulatedPotential.H"

using namespace std;

void insertionSort(int* item, double* key, int n) {
  int i, j, low;

  for (i = 1; i < n; i++) {
    low = item[i];
    j = i;
    while ((j > 0) && (key[item[j-1]] > key[low])) {
      item[j] = item[j-1];
      j--;
    }
    item[j] = low;
  }
}

double surfEnergy(double distance, double surfWidth, double barrierHeight) {
  double c = barrierHeight/(surfWidth*surfWidth);
  double x = surfWidth-distance;
  return c*x*x;
}

void writeNumbers(const char* fileName, const double* a, int n) { 
  FILE* out = fopen(fileName, "w");
    
  for (int i = 0; i < n; i++) {
    fprintf(out, "%.10g\n", a[i]);
  }
  fclose(out);
}

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
 if ( argc != 8) {
    printf("Usage: %s srcGrid potentialData.dat poreLength poreDiamX poreDiamY cornerRad outGrid\n", argv[0]);
    return 0;
  }

  const char* srcGrid = argv[1];
  const char* outGrid = argv[argc-1];

  const char* potDataFile = argv[2];

  const double poreLength = strtod(argv[3], NULL);
  const double poreDiamX = strtod(argv[4], NULL);
  const double poreDiamY = strtod(argv[5], NULL);
  const double cornerRad = strtod(argv[6], NULL);

  CylinderPore pore(poreLength, poreDiamX, poreDiamY, cornerRad);
  TabulatedPotential pot(potDataFile);

  // Load the grid.
  Grid src(srcGrid);
  printf("Loaded the grid `%s'.\n", srcGrid);

  // debug
  /*
  double v[100];
  double dr = 0.2;
  FILE* out = fopen("debug.dat", "w");
  for (int i = 0; i < 100; i++) {
    v[i] = pot.computeEnergy(i*dr);
    fprintf(out, "%.10g %.10g\n", i*dr, v[i]);
  }
  fclose(out);
  */

  double dr = pot.getInterval();
  printf("INTERVAL %g\n", dr);
  

  // Mark the pore.
  printf("Marking the pore...\n");
  const int n = src.length();
  for (int i = 0; i < n; i++) {
    Vector3 r = src.wrapDiffNearest(src.getPosition(i));
    
    if (pore.in(r)) src.setValue(i, 1.0);
    else src.setValue(i, 0.0);
  }
  src.write("test.dx");

  // Find the surface.
  printf("Marking the surface...\n");
  Grid surf(src);
  surf.blur();
  int surfNum = 0;
  for (int i = 0; i < n; i++) {
    double v = surf.getValue(i);
    if (v > 0.01 && v < 0.99) { 
      surf.setValue(i, 1.0);
      surfNum++;
    }
    else surf.setValue(i, 0.0);
  }

  // Write the surface grid.
  char outSurf[256];
  sprintf(outSurf, "surf_%s", outGrid);
  surf.write(outSurf, "surface");
  printf("Wrote the surface `%s'.\n", outSurf);

  // Make a list of surface positions.
  Vector3* surfPos = new Vector3[surfNum];
  int j = 0; 
  for (int i = 0; i < n; i++) {
    if (surf.getValue(i) > 0.0) {
      surfPos[j] = surf.getPosition(i);
      j++;
    }
  }
  printf("Extracted %d surface nodes.\n", surfNum);

  // Generate the cell decomposition.
  double cutoff = 1.01*pot.lastPosition();
  double barrierHeight = pot.computeEnergy(pot.firstPosition());
  Matrix3 sys = src.getBox();
  Vector3 origin = src.getOrigin();

  int levels = 3;
  double cutoff1 = cutoff;
  CellDecomposition** cell = new CellDecomposition*[levels];
  for (int l = 0; l < levels; l++) {
    printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], cutoff1);
    cell[l] = new CellDecomposition(sys, origin, cutoff1);
    printf("Cell decomposition of %i cells\n", cell[l]->length());
    cell[l]->decompose(surfPos, surfNum);
    //printf("Cell decomposition contains %i points\n", cell[l]->countPoints());

    cutoff1 *= 0.5;
  }
  printf("\n");

  // Now go through the grid to create the output potential.
  Grid dest(src);
  int complete = 0;
  double clockInit = omp_get_wtime();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    // If the position is inside the membrane, set it to the max value.
    if (v > 0.0) {
      dest.setValue(i, barrierHeight);
      continue;
    }

    Vector3 r = src.getPosition(i);
    IndexList neigh;
    int nNeighs = 0;

    // Find which level we need.
    // If there is anything in the neighborhood at the smallest level,
    // then we don't need to check the largest level.
    for (int l = levels-1; l >= 0; l--) {
      // Get the indices (in pos) of possible neighboring points.
       neigh = cell[l]->neighborhood(r);
       nNeighs = neigh.length();
       if (nNeighs > 0) break; // There's stuff at this level.
    }
    
    // We're at the largest level and still nothing.
    if (nNeighs < 1) {
      dest.setValue(i, 0.0);
      continue;
    }

    // Find the minimum distance to a surface point.
    Vector3 d =  src.wrapDiffNearest(surfPos[neigh.get(0)] - r);
    double minDist = d.length();
    for (int j = 1; j < nNeighs; j++) {
      Vector3 d = src.wrapDiffNearest(surfPos[neigh.get(j)] - r);
      if (d.length() < minDist) minDist = d.length();
    }
    
    // Set the value using the WallEnergy function.
    dest.setValue(i, pot.computeEnergy(minDist+2.0*dr));
    
    complete++;
    if (complete % 5000 == 0) {
      double frac = complete/double(n);
      double dur = omp_get_wtime() - clockInit;
      double remain = dur/frac - dur;
      printf("percent complete %g, time to completion %g min\n", 100.0*frac, remain/60.0);
    }
  } 

  delete[] surfPos;
  for (int l = 0; l < levels; l++) delete cell[l];
  delete[] cell;

  char comments[256];
  sprintf(comments, "potentialFile %s poreLength %g poreDiamX %g poreDiamY %g cornerRad %g", potDataFile, poreLength, poreDiamX, poreDiamY, cornerRad);
  dest.write(outGrid, comments);
  printf("Wrote `%s'.\n", outGrid);

  return 0;
}
