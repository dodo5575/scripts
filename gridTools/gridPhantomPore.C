// Author: Jeff Comer <jcomer2@illinois.edu>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "useful.H"
#include "Grid.H"
#include "CylinderPore.H"
#include "WallEnergy.H"
#include "CellDecomposition.H"

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

///////////////////////////////////////////////////////////////////////
// Driver
int main(int argc, char* argv[]) {
 if ( argc != 11 && argc != 10) {
    printf("Usage: %s srcGrid poreLength poreDiamX poreDiamY cornerRad wellWidth wellDepth barrierHeight [surfWidth] outGrid\n", argv[0]);
    printf("If `surfWidth' is specified, then the upper and lower membrane surfaces only barriers. The pore is unaffected.\n");
    return 0;
  }

  const char* srcGrid = argv[1];
  const char* outGrid = argv[argc-1];

  const double poreLength = strtod(argv[2], NULL);
  const double poreDiamX = strtod(argv[3], NULL);
  const double poreDiamY = strtod(argv[4], NULL);
  const double cornerRad = strtod(argv[5], NULL);
  const double wellWidth = strtod(argv[6], NULL);
  const double wellDepth = strtod(argv[7], NULL);
  const double barrierHeight = strtod(argv[8], NULL);
  const double surfWidth = strtod(argv[9], NULL);
  bool surfBarrier = (argc == 11);

  CylinderPore pore(poreLength, poreDiamX, poreDiamY, cornerRad);

  //double wallShift = wellWidth*1.5;
  WallEnergy wall(-fabs(wellWidth), -0.5*fabs(wellWidth), -fabs(wellDepth), 0.0, barrierHeight);
  WallEnergy wallSurf(-fabs(surfWidth), -0.5*fabs(surfWidth), 0.0, 0.0, barrierHeight);

  // Load the grid.
  Grid src(srcGrid);
  printf("Loaded the grid `%s'.\n", srcGrid);

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
  double cutoff = 1.5*wellWidth;
  Matrix3 sys = src.getBox();
  Vector3 origin = src.getOrigin();
  printf("\nGenerating cell decomposition for %s with cutoff = %.10g...\n", argv[1], cutoff);
  CellDecomposition cell(sys, origin, cutoff);
  printf("Cell decomposition of %i cells\n", cell.length());
  cell.decompose(surfPos, surfNum);
  printf("Cell decomposition contains %i points\n", cell.countPoints());

  // Now go through the grid to create the output potential.
  Grid dest(src);
  int complete = 0;
  for (int i = 0; i < n; i++) {
    double v = src.getValue(i);
    // If the position is inside the membrane, set it to the max value.
    if (v > 0.0) {
      dest.setValue(i, barrierHeight);
      continue;
    }

    // Find the distance to the surface.
    Vector3 r = src.getPosition(i);
    // Get the indices (in pos) of possible neighboring points.
    IndexList neigh = cell.neighborhood(r);
    const int nNeighs = neigh.length();
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
    if (surfBarrier) {
      double w = pore.surfaceWeight(r);
      dest.setValue(i, (1.0-w)*wall.energy(-minDist) + w*wallSurf.hardEnergy(-minDist));
    }
    else dest.setValue(i, wall.energy(-minDist));

    // Write the progress.
    int comp = (100*i)/n;
    if (abs(complete - comp) >= 5) {
      printf("%d percent complete\n", comp);
      complete = comp;
    }
  } 

  delete[] surfPos;

  char comments[256];
  sprintf(comments, "poreLength %g poreDiamX %g poreDiamY %g cornerRad %g wellWidth %g wellDepth %g barrierHeight %g", poreLength, poreDiamX, poreDiamY, cornerRad, wellWidth, wellDepth, barrierHeight);
  dest.write(outGrid, comments);
  printf("Wrote `%s'.\n", outGrid);

  return 0;
}
